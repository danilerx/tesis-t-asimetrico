# 1. LIBRERÍAS Y CONFIGURACIÓN INICIAL -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sn, dplyr, haven, progress, doSNOW, foreach, knitr, 
               ggplot2, pROC, caret, ggrepel)

# 2. CARGA Y PREPARACIÓN DE DATOS ----------------------------------------------
RUTA_CASEN <- "C:/Users/danil/Downloads/Base de datos Casen 2022 SPSS_18 marzo 2024.sav/Base de datos Casen 2022 SPSS_18 marzo 2024.sav"

columnas_estudio <- c("folio", "pco1_a", "sexo", "edad", "e6a", "o1", 
                      "ecivil", "area", "qaut", "v13", "ind_hacina", 
                      "pobreza_multi_5d", "expr", "region")

base_final <- read_sav(RUTA_CASEN, col_select = all_of(columnas_estudio)) %>%
  zap_labels() %>%
  filter(region == 13) %>%
  filter(if_all(everything(), ~ .x != -88)) %>% 
  na.omit() %>%
  group_by(folio) %>%
  summarise(
    Pobreza_Multi = as.numeric(as.character(first(pobreza_multi_5d))),
    Peso_Muestral = max(if_else(pco1_a == 1, as.numeric(expr), NA_real_), na.rm = TRUE),
    Edu_Baja = max(if_else(pco1_a == 1 & e6a <= 11, 1, 0), na.rm = TRUE),
    Desempleo = max(if_else(pco1_a == 1 & o1 != 1, 1, 0), na.rm = TRUE),
    Mujer_Jefa = max(if_else(pco1_a == 1 & sexo == 2, 1, 0), na.rm = TRUE),
    Casado_Jefe = max(if_else(pco1_a == 1 & ecivil %in% c(1, 2), 1, 0), na.rm = TRUE),
    Urbano = if_else(first(area) == 1, 1, 0),
    Bajos_Ingresos = if_else(first(qaut) <= 2, 1, 0),
    Vivienda_No_Propia = if_else(!(first(v13) %in% c(1, 2)), 1, 0),
    Hacinamiento_Ind = if_else(first(ind_hacina) > 1, 1, 0),
    .groups = 'drop'
  ) %>%
  filter(is.finite(Edu_Baja), !is.na(Peso_Muestral)) %>%
  mutate(Peso_Norm = Peso_Muestral * (n() / sum(Peso_Muestral)))

Matriz_X <- as.matrix(cbind(1, base_final %>% select(Edu_Baja:Hacinamiento_Ind)))
colnames(Matriz_X) <- c("Intercepto", "Edu_Baja", "Desempleo", "Mujer", "Casado", 
                        "Urbano", "Ingresos_Bajos", "Viv_No_Propia", "Hacinamiento")

Vector_Y <- base_final$Pobreza_Multi
Vector_Pesos <- base_final$Peso_Norm

# 3. FUNCIONES DE OPTIMIZACIÓN Y MÉTRICAS --------------------------------------
func_verosimil_t <- function(params, Y, X, pesos, nu_fijo, delta_fijo = NULL) {
  beta <- params[1:ncol(X)]
  delta <- if(is.null(delta_fijo)) params[ncol(X) + 1] else delta_fijo
  eta <- as.numeric(X %*% beta)
  prob <- pmax(pmin(pst(eta, alpha = delta, nu = nu_fijo), 1 - 1e-12), 1e-12)
  log_ind <- pesos * (Y * log(prob) + (1 - Y) * log(1 - prob))
  
  # Penalización
  v1 <- ((nu_fijo + 2) * (nu_fijo + 3)) / (3 * (nu_fijo + 1)^2)
  v2 <- 0.285 * (1 + 4 / (nu_fijo + 0.577))
  c1 <- 1 / (4 * v2); c2 <- v2 / v1
  Q_delta <- if(delta == 0) 0 else c1 * log(1 + c2 * delta^2)
  
  return(-(sum(log_ind) - Q_delta))
}

calc_aic_unificado <- function(modelo_glm, Y, X, pesos) {
  betas <- coef(modelo_glm)
  eta <- as.numeric(X %*% betas)
  prob <- switch(modelo_glm$family$link,
                 "logit" = 1 / (1 + exp(-eta)),
                 "probit" = pnorm(eta),
                 "cloglog" = 1 - exp(-exp(eta)))
  prob <- pmax(pmin(prob, 1 - 1e-12), 1e-12)
  return(-2 * sum(pesos * (Y * log(prob) + (1 - Y) * log(1 - prob))) + 2 * length(betas))
}

calc_metrics_manual <- function(fit_obj, k, n) {
  dev <- 2 * fit_obj$value
  return(c(Devianza = dev, AIC = dev + 2 * k, BIC = dev + log(n) * k))
}

calc_dev_glm <- function(mod, Y, pesos) {
  p <- pmax(pmin(predict(mod, type = "response"), 1 - 1e-12), 1e-12)
  return(-2 * sum(pesos * (Y * log(p) + (1 - Y) * log(1 - p))))
}

calc_bic_glm <- function(mod, Y, X, pesos) {
  return(calc_dev_glm(mod, Y, pesos) + log(nrow(X)) * length(coef(mod)))
}

calc_lrt <- function(ll_completo, ll_restringido, df_diff) {
  stat <- 2 * (ll_completo - ll_restringido)
  return(list(stat = stat, p = pchisq(stat, df = df_diff, lower.tail = FALSE)))
}

# 4. ESTIMACIÓN DE MODELOS Y BÚSQUEDA DE NU ÓPTIMO -----------------------------
rejilla_nu <- seq(1.44, 1.46, by = 0.001)
par_inicio <- c(rep(0, ncol(Matriz_X)), 0.1)

# Búsqueda en paralelo
cl <- makeCluster(6)
registerDoSNOW(cl)
pb <- progress_bar$new(format = " [:bar] :percent ETA: :eta", total = length(rejilla_nu), width = 60)
opts <- list(progress = function(n) pb$tick())

cat("Iniciando búsqueda de nu óptimo...\n")
lista_aic <- foreach(i = 1:length(rejilla_nu), .packages = c("sn"), 
                     .export = c("func_verosimil_t", "Vector_Y", "Matriz_X", "Vector_Pesos", "par_inicio"),
                     .options.snow = opts) %dopar% {
                       res <- optim(par = par_inicio, fn = func_verosimil_t, Y = Vector_Y, X = Matriz_X, 
                                    pesos = Vector_Pesos, nu_fijo = rejilla_nu[i], method = "BFGS")
                       return(2 * res$value + 2 * (ncol(Matriz_X) + 1))
                     }
stopCluster(cl)

nu_optimo <- rejilla_nu[which.min(unlist(lista_aic))]
cat(sprintf("Nu óptimo encontrado: %.3f\n", nu_optimo))

# Ajuste Modelos t (Asimétrico y Simétrico)
fit_t_asim <- optim(par = par_inicio, fn = func_verosimil_t, Y = Vector_Y, X = Matriz_X, 
                    pesos = Vector_Pesos, nu_fijo = nu_optimo, method = "BFGS", hessian = TRUE)

fit_t_sim <- optim(par = rep(0, ncol(Matriz_X)), fn = func_verosimil_t, Y = Vector_Y, X = Matriz_X, 
                   pesos = Vector_Pesos, nu_fijo = nu_optimo, delta_fijo = 0, method = "BFGS")

# Ajuste Modelos Clásicos GLM
mod_logit <- glm(Vector_Y ~ Matriz_X[,-1], family = binomial(link = "logit"), weights = Vector_Pesos)
mod_probit <- glm(Vector_Y ~ Matriz_X[,-1], family = binomial(link = "probit"), weights = Vector_Pesos)
mod_cloglog <- glm(Vector_Y ~ Matriz_X[,-1], family = binomial(link = "cloglog"), weights = Vector_Pesos)

# 5. DIAGNÓSTICO Y TABLAS DE RESULTADOS ----------------------------------------
n <- nrow(Matriz_X)
m_tasim <- calc_metrics_manual(fit_t_asim, length(fit_t_asim$par), n)
m_tsim <- calc_metrics_manual(fit_t_sim, length(fit_t_sim$par), n)

tabla_comparativa <- data.frame(
  Modelo = c("t-Asimétrico", "t-Simétrico", "Logit", "Probit", "Cloglog"),
  Devianza = c(m_tasim["Devianza"], m_tsim["Devianza"], calc_dev_glm(mod_logit, Vector_Y, Vector_Pesos),
               calc_dev_glm(mod_probit, Vector_Y, Vector_Pesos), calc_dev_glm(mod_cloglog, Vector_Y, Vector_Pesos)),
  AIC = c(m_tasim["AIC"], m_tsim["AIC"], calc_aic_unificado(mod_logit, Vector_Y, Matriz_X, Vector_Pesos),
          calc_aic_unificado(mod_probit, Vector_Y, Matriz_X, Vector_Pesos), calc_aic_unificado(mod_cloglog, Vector_Y, Matriz_X, Vector_Pesos)),
  BIC = c(m_tasim["BIC"], m_tsim["BIC"], calc_bic_glm(mod_logit, Vector_Y, Matriz_X, Vector_Pesos),
          calc_bic_glm(mod_probit, Vector_Y, Matriz_X, Vector_Pesos), calc_bic_glm(mod_cloglog, Vector_Y, Matriz_X, Vector_Pesos))
) %>% arrange(AIC)

# Errores estándar y p-valores del modelo propuesto
ee <- sqrt(diag(solve(fit_t_asim$hessian)))
reporte_final <- data.frame(
  Variable = c(colnames(Matriz_X), "Delta"),
  Estimacion = round(fit_t_asim$par, 4),
  Error_Std = round(ee, 4),
  P_Valor = round(2 * (1 - pnorm(abs(fit_t_asim$par / ee))), 6)
)

print(kable(tabla_comparativa, digits = 2, caption = "Comparación de Criterios de Información"))
print(kable(reporte_final, caption = "Estimaciones del Modelo t-Asimétrico"))

# Prueba LRT para asimetría
lrt_asimetria <- calc_lrt(ll_completo = -fit_t_asim$value, ll_restringido = -fit_t_sim$value, df_diff = 1)
cat(sprintf("\nPrueba LRT (Asimétrico vs Simétrico): Estadístico = %.3f, p-valor = %.6f\n", lrt_asimetria$stat, lrt_asimetria$p))

# 6. GRÁFICOS PRINCIPALES DE LA TESIS ------------------------------------------

# 6.1 Curva ROC Formal
betas_opt <- fit_t_asim$par[1:ncol(Matriz_X)]
delta_opt <- fit_t_asim$par[ncol(Matriz_X) + 1]
eta_tasim <- as.numeric(Matriz_X %*% betas_opt)
p_hat <- pmax(pmin(pst(eta_tasim, alpha = delta_opt, nu = nu_optimo), 1 - 1e-12), 1e-12)

roc_tasim <- roc(Vector_Y, p_hat, quiet = TRUE)
roc_logit <- roc(Vector_Y, predict(mod_logit, type = "response"), quiet = TRUE)

df_roc <- rbind(
  data.frame(especificidad = roc_tasim$specificities, sensibilidad = roc_tasim$sensitivities, modelo = "t-Asimétrico"),
  data.frame(especificidad = roc_logit$specificities, sensibilidad = roc_logit$sensitivities, modelo = "Logit")
)

p_roc <- ggplot(df_roc, aes(x = especificidad, y = sensibilidad, color = modelo, linetype = modelo)) +
  scale_x_reverse(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey80", size = 0.5) +
  geom_line(size = 0.8) +
  scale_color_manual(values = c("t-Asimétrico" = "#001f3f", "Logit" = "#4F4F4F")) +
  scale_linetype_manual(values = c("t-Asimétrico" = "solid", "Logit" = "dashed")) +
  labs(title = "Comparación de Curvas ROC",
       subtitle = sprintf("AUC t-Asimétrico: %.3f | AUC Logit: %.3f", auc(roc_tasim), auc(roc_logit)),
       x = "Especificidad", y = "Sensibilidad") +
  theme_classic() + 
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

ggsave("ROC_Comparativa.pdf", p_roc, width = 6, height = 5)

# 6.2 Visualización del Trade-off (Verosimilitud Penalizada)
comparar_penalizacion <- function(d, Y, X, pesos, betas, nu) {
  eta <- as.numeric(X %*% betas)
  prob <- pmax(pmin(sn::pst(eta, alpha = d, nu = nu), 1 - 1e-12), 1e-12)
  ll_ind <- sum(pesos * (Y * log(prob) + (1 - Y) * log(1 - prob)))
  
  v1 <- ((nu + 2) * (nu + 3)) / (3 * (nu + 1)^2)
  v2 <- 0.285 * (1 + 4 / (nu + 0.577))
  Q_delta <- if(d == 0) 0 else (1 / (4 * v2)) * log(1 + (v2 / v1) * d^2)
  
  return(data.frame(Delta = d, PLL = ll_ind - Q_delta))
}

rango_delta <- seq(0, 5, length.out = 100)
df_curvas <- purrr::map_df(rango_delta, ~comparar_penalizacion(.x, Vector_Y, Matriz_X, Vector_Pesos, betas_opt, nu_optimo))
delta_pll <- df_curvas$Delta[which.max(df_curvas$PLL)]

p_penalizacion <- ggplot(df_curvas, aes(x = Delta, y = PLL)) +
  geom_line(size = 1.1, color = "#001f3f") +
  geom_vline(xintercept = delta_pll, linetype = "dashed", color = "#e74c3c") +
  annotate("label", x = delta_pll, y = max(df_curvas$PLL), 
           label = paste0("hat(delta) == ", round(delta_pll, 4)), 
           parse = TRUE, color = "white", fill = "#e74c3c", fontface = "bold") +
  labs(title = "Log-Verosimilitud Penalizada Perfilada",
       x = expression(delta), y = "Log-Verosimilitud Penalizada") +
  theme_minimal()

ggsave("Verosimilitud_Penalizada.pdf", p_penalizacion, width = 6, height = 5)