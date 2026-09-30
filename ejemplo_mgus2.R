# =============================================================================
# Ejemplo de 5 minutos: progresión de MGUS a neoplasia de células plasmáticas
# (PCM) con la muerte sin PCM como riesgo competitivo.
#
# Datos públicos: survival::mgus2, 1384 pacientes de la Mayo Clinic
# (Kyle et al., NEJM 2002). Viene con el paquete survival.
#
# Ejecuta desde la carpeta del repositorio:   Rscript ejemplo_mgus2.R
# Resultados en resultados/: figura PNG, tabla AJ frente a 1 - KM, modelos y
# metadatos. Uso docente; no apto para decisiones clínicas.
# =============================================================================

# Trabaja en la carpeta del script aunque se llame desde otra (Rscript ruta/ejemplo_mgus2.R).
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
# En Linux y Mac, Rscript codifica como "~+~" los espacios de la ruta.
script <- gsub("~+~", " ", script, fixed = TRUE)
if (length(script) == 1) setwd(dirname(normalizePath(script)))
if (!file.exists("R/riesgos_competitivos.R")) {
  stop("No encuentro R/riesgos_competitivos.R. Ejecuta el script desde la carpeta del repositorio.")
}
source("R/riesgos_competitivos.R")

# 1. Datos: construcción estándar de la viñeta "compete" de survival -----------
mgus2 <- survival::mgus2
mgus2$etime <- ifelse(mgus2$pstat == 0, mgus2$futime, mgus2$ptime)  # meses
mgus2$event <- ifelse(mgus2$pstat == 0, 2 * mgus2$death, 1)         # 0 censura, 1 PCM, 2 muerte
mgus2$years <- mgus2$etime / 12  # en años se lee mejor; no cambia ninguna estimación

datos <- prepare_competing_data(
  mgus2, time = "years", status = "event",
  event_codes = 1, competing_codes = 2, censor_codes = 0,
  covariates = c("age", "sex"),
  event_label = "PCM", competing_label = "muerte sin PCM"
)

# 2. Parámetros ----------------------------------------------------------------
times <- c(5, 10, 15, 20, 25, 30)  # años
min_at_risk <- 20                  # por debajo, la estimación se marca como inestable
covariates <- c("age", "sex")

# 3. Análisis ------------------------------------------------------------------
tabla <- compare_aj_km(datos, times, min_at_risk = min_at_risk)
modelos <- fit_competing_models(datos, covariates)

plot_aj_vs_km(
  datos, file = "resultados/mgus2_figura.png", time_unit = "años",
  xmax = 30, xticks = seq(0, 30, 5), highlight_time = 20, min_at_risk = min_at_risk,
  main = "Progresión de MGUS a neoplasia de células plasmáticas (PCM)",
  caption = paste("Datos: survival::mgus2 (Kyle et al., NEJM 2002).",
                  "Uso docente; no apto para decisiones clínicas.")
)
write_table(tabla, "resultados/mgus2_tabla_aj_km.csv")
write_table(modelos, "resultados/mgus2_modelos.csv")
write_metadata(
  "resultados/mgus2_metadatos.txt", datos,
  input = paste0("survival::mgus2 (survival ", utils::packageDescription("survival")$Version,
                 "); tiempo = meses / 12; estado: 0 censura, 1 PCM, 2 muerte sin PCM"),
  parameters = list(unidad_tiempo = "años", tiempos = times, min_en_riesgo = min_at_risk,
                    covariables = covariates, script = "ejemplo_mgus2.R")
)

# 4. Resumen -------------------------------------------------------------------
cat("\nMGUS (survival::mgus2): n = ", nrow(datos), "; ", sum(datos$status == 1), " PCM, ",
    sum(datos$status == 2), " muertes sin PCM y ", sum(datos$status == 0), " censuras.\n\n",
    "Incidencia acumulada de PCM (tiempo en años): ",
    "Aalen-Johansen (AJ) frente a 1 - Kaplan-Meier (1-KM)\n",
    sep = "")
cat(format_comparison(tabla), sep = "\n")
at20 <- tabla[tabla$time == 20, ]
cat("\nA 20 años, 1-KM da ", format_num(100 * at20$one_minus_km, 1), " % y AJ ",
    format_num(100 * at20$aj, 1), " %: 1-KM sobreestima ", format_num(at20$diff_pp, 1),
    " puntos (x", format_num(at20$ratio, 1), ").\n\n", sep = "")
cat("Modelos con edad y sexo (HR e IC 95 %; n = ", modelos$n[1], "):\n", sep = "")
cat(format_models(modelos), sep = "\n")
cat("\nFicheros en resultados/: mgus2_figura.png, mgus2_tabla_aj_km.csv,",
    "mgus2_modelos.csv y mgus2_metadatos.txt\n")
