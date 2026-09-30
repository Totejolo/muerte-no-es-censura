# =============================================================================
# Plantilla: incidencia acumulada con riesgos competitivos para tus datos (CSV).
#
# 1. Copia tu CSV en la carpeta datos/ (git la ignora: tus datos no se suben).
# 2. Edita solo el bloque CONFIGURACIÓN.
# 3. Ejecuta desde la carpeta del repositorio:   Rscript plantilla.R
#
# Formato del CSV: una fila por paciente, con el tiempo desde el origen (p. ej.,
# el diagnóstico) hasta el primer evento o la censura, y el estado en ese momento.
# Sin editar nada, funciona con un CSV de ejemplo que se genera desde
# survival::mgus2 (datos públicos) la primera vez que se ejecuta.
# Uso docente y de investigación; no apto para decisiones clínicas.
# =============================================================================

# ---- CONFIGURACIÓN ----------------------------------------------------------
input_csv       <- "datos/ejemplo_mgus2.csv"  # tu CSV (ruta relativa a esta carpeta o absoluta)
csv_sep         <- ","             # separador de columnas; Excel en español suele usar ";"
csv_dec         <- "."             # separador decimal; con ";" suele ser ","
csv_encoding    <- "UTF-8-BOM"     # lee UTF-8 con o sin BOM; CSV antiguos de Excel: "latin1"
time_col        <- "meses"         # columna con el tiempo de seguimiento
status_col      <- "estado"        # columna con el estado al final del seguimiento
censor_codes    <- "censura"       # códigos de censura; pueden ser varios: c("0", "perdido")
event_codes     <- "PCM"           # códigos del evento de interés
competing_codes <- "muerte"        # códigos de los eventos competidores
event_label     <- "PCM"           # nombres cortos para el gráfico y las tablas
competing_label <- "muerte sin PCM"
time_unit       <- "meses"         # unidad del tiempo (ejes y metadatos)
times           <- c(60, 120, 180, 240)  # tiempos de la tabla, en esa unidad
covariates      <- c("edad", "sexo")     # covariables de los modelos; character(0) = sin modelos
min_at_risk     <- 20              # por debajo, la estimación se marca como inestable
plot_xmax       <- NULL            # final del eje X; NULL = todo el seguimiento
output_prefix   <- "resultados/plantilla"  # prefijo de los ficheros de salida
# -----------------------------------------------------------------------------

# Trabaja en la carpeta del script aunque se llame desde otra (Rscript ruta/plantilla.R).
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
# En Linux y Mac, Rscript codifica como "~+~" los espacios de la ruta.
script <- gsub("~+~", " ", script, fixed = TRUE)
if (length(script) == 1) setwd(dirname(normalizePath(script)))
if (!file.exists("R/riesgos_competitivos.R")) {
  stop("No encuentro R/riesgos_competitivos.R. Ejecuta el script desde la carpeta del repositorio.")
}
source("R/riesgos_competitivos.R")

# CSV de ejemplo: se genera en cada equipo desde survival::mgus2, así el
# repositorio no guarda copias de datos de terceros.
example_csv <- "datos/ejemplo_mgus2.csv"
if (identical(input_csv, example_csv) && !file.exists(example_csv)) {
  m <- survival::mgus2
  example <- data.frame(
    meses = ifelse(m$pstat == 0, m$futime, m$ptime),
    estado = ifelse(m$pstat == 1, "PCM", ifelse(m$death == 1, "muerte", "censura")),
    edad = m$age,
    sexo = m$sex
  )
  dir.create("datos", showWarnings = FALSE)
  write.csv(example, example_csv, row.names = FALSE)
  message("Creado ", example_csv, " desde survival::mgus2 (datos públicos de ejemplo).")
}
if (!file.exists(input_csv)) stop("No existe el fichero ", input_csv, ". Revisa 'input_csv'.")

# 1. Leer y validar ------------------------------------------------------------
# Un CSV en latin1 leído como UTF-8 se corta a medias y el error que sale después
# despista; se comprueba antes.
if (grepl("^UTF-?8", toupper(csv_encoding)) &&
    !all(validUTF8(readLines(input_csv, warn = FALSE)))) {
  stop("El fichero ", input_csv, " no está en UTF-8. Prueba con csv_encoding <- \"latin1\" ",
       "(CSV antiguos de Excel).")
}
raw <- utils::read.csv(input_csv, sep = csv_sep, dec = csv_dec, fileEncoding = csv_encoding,
                       stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("", "NA"))
if (nrow(raw) == 0 || ncol(raw) < 2) {
  stop("No se ha podido leer bien ", input_csv, " (", nrow(raw), " filas, ", ncol(raw),
       " columnas). Revisa csv_sep, csv_dec y csv_encoding.")
}
datos <- prepare_competing_data(
  raw, time = time_col, status = status_col,
  event_codes = event_codes, competing_codes = competing_codes, censor_codes = censor_codes,
  covariates = covariates, event_label = event_label, competing_label = competing_label
)

# 2. Tabla, figura y modelos ---------------------------------------------------
tabla <- compare_aj_km(datos, times, min_at_risk = min_at_risk)
# En la figura se destaca el último tiempo de la tabla con suficientes pacientes en riesgo.
reliable <- tabla$time[!tabla$few_at_risk]
highlight <- if (length(reliable) > 0) max(reliable) else NULL

files <- paste0(output_prefix,
                c("_figura.png", "_tabla_aj_km.csv", "_modelos.csv", "_metadatos.txt"))
plot_aj_vs_km(datos, file = files[1], time_unit = time_unit, xmax = plot_xmax,
              highlight_time = highlight, min_at_risk = min_at_risk,
              caption = paste0("Datos: ", basename(input_csv),
                               ". Uso docente; no apto para decisiones clínicas."))
write_table(tabla, files[2], sep = csv_sep, dec = csv_dec)
modelos <- NULL
if (length(covariates) > 0) {
  modelos <- fit_competing_models(datos, covariates)
  write_table(modelos, files[3], sep = csv_sep, dec = csv_dec)
}
write_metadata(
  files[4], datos,
  input = paste0(input_csv, " (md5 ", unname(tools::md5sum(input_csv)), ")"),
  parameters = list(
    columna_tiempo = time_col, columna_estado = status_col,
    codigos_censura = censor_codes, codigos_evento = event_codes,
    codigos_competidor = competing_codes, unidad_tiempo = time_unit, tiempos = times,
    covariables = if (length(covariates) > 0) covariates else "ninguna",
    min_en_riesgo = min_at_risk, script = "plantilla.R"
  )
)

# 3. Resumen -------------------------------------------------------------------
cat("\n", basename(input_csv), ": n = ", nrow(datos), " (", event_label, ": ",
    sum(datos$status == 1), "; ", competing_label, ": ", sum(datos$status == 2),
    "; censura: ", sum(datos$status == 0), ")\n\n", "Incidencia acumulada de ", event_label,
    " (tiempo en ", time_unit, "): Aalen-Johansen (AJ) frente a 1 - Kaplan-Meier (1-KM)\n",
    sep = "")
cat(format_comparison(tabla), sep = "\n")
if (!is.null(modelos)) {
  cat("\nModelos (HR e IC 95 %; n = ", modelos$n[1], "):\n", sep = "")
  cat(format_models(modelos), sep = "\n")
}
written <- files[c(TRUE, TRUE, !is.null(modelos), TRUE)]
cat("\nFicheros:", paste(written, collapse = ", "), "\n")
