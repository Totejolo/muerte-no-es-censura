# =============================================================================
# riesgos_competitivos.R
#
# Incidencia acumulada con riesgos competitivos: Aalen-Johansen frente a
# 1 - Kaplan-Meier, Fine-Gray (sHR) y Cox por causa (csHR).
# Solo usa base R y survival (paquete recomendado que viene con R).
# Uso docente y de investigación; no sirve para decisiones clínicas.
#
# Funciones principales
#   prepare_competing_data()  valida los datos y recodifica el estado a 0/1/2
#   compare_aj_km()           tabla AJ (IC 95 %) frente a 1 - KM a tiempos dados
#   plot_aj_vs_km()           gráfico AJ frente a 1 - KM, en pantalla o en PNG
#   fit_competing_models()    Fine-Gray y Cox por causa en una tabla ordenada
# Auxiliares
#   format_comparison(), format_models()  resúmenes de texto en español
#   write_table(), write_metadata()       exportación con metadatos
# =============================================================================

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Falta el paquete 'survival'. Instálalo con install.packages(\"survival\").")
}
library(survival)

rc_version <- "1.0.0"

# Tras prepare_competing_data() el estado queda codificado así:
#   0 = censura, 1 = evento de interés, 2 = evento competidor.
# Estos nombres se usan por dentro y no valen como covariables.
reserved_names <- c("time", "status", "subject_id",
                    "fgstart", "fgstop", "fgstatus", "fgwt")


# -----------------------------------------------------------------------------
# Preparación y validación
# -----------------------------------------------------------------------------

#' Valida los datos y devuelve un data.frame con `time`, `status` (0/1/2) y
#' las covariables pedidas. Los códigos pueden ser numéricos o de texto y cada
#' lista admite varios códigos (p. ej., competing_codes = c("muerte", "trasplante")).
#' Las etiquetas se usan en el gráfico y en la tabla de modelos.
prepare_competing_data <- function(data, time, status,
                                   event_codes, competing_codes, censor_codes,
                                   covariates = character(0),
                                   event_label = "evento",
                                   competing_label = "competidor") {
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' debe ser un data.frame con al menos una fila.", call. = FALSE)
  }
  absent <- setdiff(c(time, status, covariates), names(data))
  if (length(absent) > 0) {
    stop("No existen estas columnas: ", paste(absent, collapse = ", "),
         ". Columnas disponibles: ", paste(names(data), collapse = ", "), ".",
         call. = FALSE)
  }
  clash <- intersect(covariates, reserved_names)
  if (length(clash) > 0) {
    stop("Renombra estas covariables (el nombre se usa por dentro): ",
         paste(clash, collapse = ", "), ".", call. = FALSE)
  }
  invalid <- covariates[make.names(covariates) != covariates]
  if (length(invalid) > 0) {
    stop("Nombres de covariable no válidos en R: ", paste(invalid, collapse = ", "),
         ". Usa letras, números, '.' o '_', sin espacios.", call. = FALSE)
  }

  # Tiempo: numérico, sin NA, finito y no negativo.
  tt <- data[[time]]
  if (!is.numeric(tt)) {
    stop("La columna de tiempo '", time, "' no es numérica. Si viene de un CSV, ",
         "revisa el separador decimal.", call. = FALSE)
  }
  if (anyNA(tt)) stop(rows_msg("Hay NA", time, is.na(tt)), call. = FALSE)
  if (any(is.infinite(tt))) {
    stop(rows_msg("Hay tiempos infinitos", time, is.infinite(tt)), call. = FALSE)
  }
  if (any(tt < 0)) stop(rows_msg("Hay tiempos negativos", time, tt < 0), call. = FALSE)

  # Estado: se compara como texto sin espacios, así 1, "1" y " 1" son el mismo código.
  st <- data[[status]]
  if (anyNA(st)) stop(rows_msg("Hay NA", status, is.na(st)), call. = FALSE)
  st <- trimws(as.character(st))
  codes <- list(censura = censor_codes, evento = event_codes, competidor = competing_codes)
  codes <- lapply(codes, function(x) unique(trimws(as.character(x))))
  for (k in names(codes)) {
    if (length(codes[[k]]) == 0 || anyNA(codes[[k]])) {
      stop("Declara al menos un código de ", k, " (sin NA).", call. = FALSE)
    }
  }
  all_codes <- unlist(codes, use.names = FALSE)
  repeated <- unique(all_codes[duplicated(all_codes)])
  if (length(repeated) > 0) {
    stop("Códigos repetidos en más de una lista (censura, evento, competidor): ",
         paste(repeated, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(unique(st), all_codes)
  if (length(unknown) > 0) {
    counts <- table(st[st %in% unknown])
    stop("Códigos de estado no declarados en '", status, "': ",
         paste0("'", names(counts), "' (", as.integer(counts), " filas)", collapse = ", "),
         ". Decláralos como censura, evento o competidor.", call. = FALSE)
  }

  new_status <- ifelse(st %in% codes$evento, 1L, ifelse(st %in% codes$competidor, 2L, 0L))
  if (!any(new_status == 1L)) {
    stop("No hay ningún evento de interés (códigos: ",
         paste(codes$evento, collapse = ", "), "). Revisa los códigos.", call. = FALSE)
  }
  if (!any(new_status == 2L)) {
    message("No hay eventos competidores: Aalen-Johansen coincidirá con 1 - Kaplan-Meier.")
  }

  out <- data.frame(time = tt, status = new_status)
  for (v in covariates) {
    x <- data[[v]]
    # Los NA se admiten (los modelos usan casos completos); los infinitos, no.
    if (is.numeric(x) && any(is.infinite(x))) {
      stop(rows_msg("Hay valores infinitos", v, is.infinite(x)), call. = FALSE)
    }
    out[[v]] <- x
  }
  attr(out, "event_label") <- event_label
  attr(out, "competing_label") <- competing_label
  out
}

# Mensaje de error que señala las primeras filas problemáticas.
rows_msg <- function(what, column, bad) {
  rows <- which(bad)
  shown <- paste(utils::head(rows, 5), collapse = ", ")
  if (length(rows) > 5) shown <- paste0(shown, " y ", length(rows) - 5, " más")
  paste0(what, " en '", column, "' (filas: ", shown, ").")
}

check_min_at_risk <- function(min_at_risk) {
  ok <- is.numeric(min_at_risk) && length(min_at_risk) == 1 && is.finite(min_at_risk) &&
    min_at_risk >= 0
  if (!ok) stop("'min_at_risk' debe ser un único número finito >= 0.", call. = FALSE)
}

check_prepared <- function(prepared) {
  ok <- is.data.frame(prepared) && all(c("time", "status") %in% names(prepared)) &&
    is.numeric(prepared$time) && all(prepared$status %in% 0:2)
  if (!ok) stop("Prepara antes los datos con prepare_competing_data().", call. = FALSE)
}

# Etiqueta guardada por prepare_competing_data(). Al subdividir un data.frame se
# pierden sus atributos; entonces se usan las etiquetas por defecto.
label_of <- function(prepared, which) {
  x <- attr(prepared, paste0(which, "_label"))
  if (is.null(x)) c(event = "evento", competing = "competidor")[[which]] else x
}

# Curvas de las que sale todo lo demás:
#   aj: Aalen-Johansen (survfit con estado multinivel);
#   km: Kaplan-Meier del evento de interés tratando al competidor como censura.
fit_curves <- function(prepared) {
  state <- factor(prepared$status, levels = 0:2, labels = c("censored", "event", "competing"))
  aj <- survfit(Surv(prepared$time, state) ~ 1)
  km <- survfit(Surv(prepared$time, prepared$status == 1) ~ 1)
  list(aj = aj, km = km, j = match("event", aj$states))
}

# Pacientes en riesgo en t: los que siguen en seguimiento (time >= t).
# Equivale a sum(time >= t) para cada t, pero sin recorrer los datos una vez por tiempo.
n_at_risk <- function(prepared, times) {
  length(prepared$time) - findInterval(times, sort(prepared$time), left.open = TRUE)
}


# -----------------------------------------------------------------------------
# Tabla comparativa
# -----------------------------------------------------------------------------

#' Incidencia acumulada del evento de interés a los tiempos pedidos:
#' en riesgo, AJ con IC 95 % (el de survfit), 1 - KM, diferencia en puntos
#' porcentuales, razón (1 - KM) / AJ y aviso de pocos en riesgo.
#' Más allá del último seguimiento (nadie en riesgo) devuelve NA.
compare_aj_km <- function(prepared, times, min_at_risk = 20) {
  check_prepared(prepared)
  check_min_at_risk(min_at_risk)
  if (!is.numeric(times) || length(times) == 0 || any(!is.finite(times)) || any(times < 0)) {
    stop("'times' debe ser un vector numérico de tiempos finitos >= 0, sin NA.", call. = FALSE)
  }
  times <- sort(unique(as.numeric(times)))
  curves <- fit_curves(prepared)
  aj <- summary(curves$aj, times = times, extend = TRUE)
  km <- summary(curves$km, times = times, extend = TRUE)
  # summary() devuelve una fila por tiempo y una columna por estado; si la forma no
  # cuadra, mejor parar que reciclar valores de otro estado.
  column <- function(m) {
    if (length(m) != length(times) * length(curves$aj$states)) {
      stop("survfit no ha devuelto una estimación por cada tiempo pedido.", call. = FALSE)
    }
    matrix(m, nrow = length(times))[, curves$j]
  }

  tab <- data.frame(
    time = times,
    n_risk = n_at_risk(prepared, times),
    aj = column(aj$pstate),
    aj_lower = column(aj$lower),
    aj_upper = column(aj$upper),
    one_minus_km = 1 - km$surv
  )
  # Sin nadie en seguimiento no hay estimación: NA en vez de arrastrar el último valor.
  beyond <- tab$n_risk == 0
  tab[beyond, c("aj", "aj_lower", "aj_upper", "one_minus_km")] <- NA
  tab$diff_pp <- 100 * (tab$one_minus_km - tab$aj)
  tab$ratio <- ifelse(!is.na(tab$aj) & tab$aj > 0, tab$one_minus_km / tab$aj, NA)
  tab$few_at_risk <- tab$n_risk < min_at_risk

  few <- tab$time[tab$few_at_risk & !beyond]
  if (length(few) > 0) {
    message("Pocos pacientes en riesgo (< ", min_at_risk, ") en t = ",
            paste(format_time(few), collapse = ", "), ": estimaciones inestables.")
  }
  if (any(beyond)) {
    message("Nadie sigue en seguimiento en t = ",
            paste(format_time(tab$time[beyond]), collapse = ", "),
            ": se devuelve NA.")
  }
  tab
}


# -----------------------------------------------------------------------------
# Gráfico
# -----------------------------------------------------------------------------

#' Dibuja AJ frente a 1 - KM (eje en %), con la tabla de pacientes en riesgo bajo
#' el eje y una zona gris donde quedan menos de `min_at_risk`. Con `file` escribe
#' un PNG; sin `file` dibuja en el dispositivo activo. `highlight_time` marca la
#' diferencia entre curvas en ese tiempo.
plot_aj_vs_km <- function(prepared, file = NULL, time_unit = "meses",
                          xmax = NULL, xticks = NULL, highlight_time = NULL,
                          min_at_risk = 20, main = NULL, caption = NULL,
                          width = 2000, height = 1350, res = 250) {
  check_prepared(prepared)
  check_min_at_risk(min_at_risk)
  if (is.null(xmax)) xmax <- max(prepared$time)
  if (!is.numeric(xmax) || length(xmax) != 1 || !is.finite(xmax) || xmax <= 0) {
    stop("'xmax' debe ser un número finito mayor que 0.", call. = FALSE)
  }
  xy <- curve_points(prepared, xmax)
  aj_xy <- xy$aj
  km_xy <- xy$km
  if (is.null(xticks)) xticks <- pretty(c(0, xmax))
  xticks <- xticks[xticks >= 0 & xticks <= xmax]
  yticks <- pretty(c(0, km_xy$y, aj_xy$y))
  if (is.null(main)) main <- paste0("Incidencia acumulada: ", label_of(prepared, "event"))

  # Colores Okabe-Ito, aptos para daltonismo; el tipo de línea es una segunda pista.
  col_aj <- "#0072B2"
  col_km <- "#D55E00"
  ink <- "#1F1F1F"
  ink2 <- "#52514E"

  # Con locale C (p. ej., en contenedores) cairo no dibuja tildes: se usa UTF-8
  # solo mientras se dibuja.
  restore_ctype <- use_utf8_ctype()
  on.exit(restore_ctype(), add = TRUE)
  if (!is.null(file)) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    png_args <- list(filename = file, width = width, height = height, res = res)
    # En macOS, quartz: cairo necesita XQuartz y, si falta, png() solo avisa y no
    # escribe nada. En Linux y Windows, cairo da el mismo dibujo en los dos.
    if (isTRUE(capabilities("aqua"))) {
      png_args$type <- "quartz"
    } else if (isTRUE(capabilities("cairo"))) {
      png_args$type <- "cairo"
    }
    before <- grDevices::dev.cur()
    do.call(grDevices::png, png_args)
    if (identical(grDevices::dev.cur(), before)) {
      stop("No se ha podido abrir el dispositivo PNG para escribir ", file,
           ". Revisa el soporte gráfico de R con capabilities().", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  op <- graphics::par(mar = c(7.6, 6.2, 4.4, 1.6), mgp = c(3, 0.6, 0), tcl = -0.25,
                      las = 1, family = "sans", col.axis = ink2)
  on.exit(graphics::par(op), add = TRUE, after = FALSE)

  plot(NA, xlim = c(0, xmax), ylim = range(yticks), yaxs = "i", axes = FALSE, ann = FALSE)
  usr <- graphics::par("usr")

  legend_args <- list("topleft", inset = c(0.01, 0.01), bg = "white", box.col = NA, cex = 0.76,
                      text.col = ink, seg.len = 2.6, y.intersp = 1.25,
                      legend = c("1 − Kaplan-Meier: trata el riesgo competitivo como censura",
                                 "Aalen-Johansen: tiene en cuenta el riesgo competitivo"),
                      col = c(col_km, col_aj), lty = c("42", "solid"), lwd = c(2, 2.2))
  legend_box <- do.call(graphics::legend, c(legend_args, plot = FALSE))$rect

  # Zona gris: menos de min_at_risk pacientes en riesgo (cola inestable).
  ut <- sort(unique(prepared$time))
  k <- which(n_at_risk(prepared, ut) < min_at_risk)[1]
  if (!is.na(k)) {
    x0 <- if (k > 1) ut[k - 1] else 0
    if (x0 < xmax) {
      graphics::rect(x0, usr[3], usr[2], usr[4], col = "#EFEEEA", border = NA)
      # Etiqueta dentro de la zona; si no cabe, alineada al borde derecho. Si la
      # zona empieza bajo la leyenda, la etiqueta va debajo para que no la tape.
      label <- paste0("pocos en riesgo\n(< ", min_at_risk, ")")
      fits <- x0 + 1.1 * graphics::strwidth(label, cex = 0.62) <= usr[2]
      under_legend <- x0 < legend_box$left + legend_box$w
      graphics::text(if (fits) x0 else usr[2],
                     if (under_legend) legend_box$top - legend_box$h else usr[4], label,
                     adj = c(if (fits) -0.05 else 1.05, 1.25), cex = 0.62, col = ink2)
    }
  }
  graphics::abline(h = yticks, col = "#E4E3DD", lwd = 0.8)
  graphics::abline(h = 0, col = "#B9B8B0", lwd = 1)
  graphics::axis(2, at = yticks, labels = paste0(prettyNum(yticks, decimal.mark = ","), " %"),
                 lwd = 0, cex.axis = 0.8)
  graphics::axis(1, at = xticks, labels = prettyNum(xticks, decimal.mark = ","),
                 lwd = 0, lwd.ticks = 1, col.ticks = "#B9B8B0", cex.axis = 0.8)

  graphics::lines(km_xy$x, km_xy$y, type = "s", col = col_km, lwd = 2, lty = "42")
  graphics::lines(aj_xy$x, aj_xy$y, type = "s", col = col_aj, lwd = 2.2)

  if (!is.null(highlight_time) && highlight_time <= max(aj_xy$x)) {
    v_aj <- aj_xy$y[findInterval(highlight_time, aj_xy$x)]
    v_km <- km_xy$y[findInterval(highlight_time, km_xy$x)]
    graphics::segments(highlight_time, v_aj, highlight_time, v_km, col = ink2, lwd = 1)
    graphics::points(rep(highlight_time, 2), c(v_aj, v_km), pch = 21, cex = 1.3, lwd = 1.5,
                     col = "white", bg = c(col_aj, col_km))
    # Valores centrados: 1 - KM encima de su punto y AJ debajo del suyo.
    graphics::text(highlight_time, v_km, paste0(format_num(v_km, 1), " %"),
                   adj = c(0.5, -0.9), cex = 0.72, col = ink)
    graphics::text(highlight_time, v_aj, paste0(format_num(v_aj, 1), " %"),
                   adj = c(0.5, 1.9), cex = 0.72, col = ink)
    ratio <- if (v_aj > 0) paste0("\n(×", format_num(v_km / v_aj, 1), ")") else ""
    graphics::text(highlight_time, (v_aj + v_km) / 2,
                   paste0("+", format_num(v_km - v_aj, 1), " pp", ratio),
                   pos = 4, cex = 0.72, col = ink)
  }

  do.call(graphics::legend, legend_args)

  graphics::title(main = main, adj = 0, line = 2.3, cex.main = 1, col.main = ink)
  graphics::mtext(paste0("Riesgo competitivo: ", label_of(prepared, "competing"),
                         "  ·  n = ", nrow(prepared)),
                  side = 3, line = 0.9, adj = 0, cex = 0.78, col = ink2)
  graphics::mtext("Incidencia acumulada", side = 2, line = 3.4, las = 0, cex = 0.82, col = ink)
  graphics::mtext(paste0("Tiempo (", time_unit, ")"), side = 1, line = 1.9, cex = 0.82, col = ink)
  graphics::mtext("En riesgo", side = 1, line = 3.7, at = usr[1], adj = 1, cex = 0.72,
                  col = ink, font = 2)
  graphics::mtext(n_at_risk(prepared, xticks), side = 1, line = 3.7, at = xticks,
                  cex = 0.72, col = ink2)
  if (!is.null(caption)) {
    graphics::mtext(caption, side = 1, line = 5.9, at = usr[1], adj = 0, cex = 0.62, col = ink2)
  }
  invisible(file)
}

# Puntos de las dos curvas escalonadas, en %. Se cortan en xmax para no dibujar
# fuera del eje y en el último seguimiento para no prolongar la curva donde
# compare_aj_km() devuelve NA.
curve_points <- function(prepared, xmax) {
  curves <- fit_curves(prepared)
  end <- min(xmax, max(prepared$time))
  cut_at <- function(x, y) {
    keep <- x <= end
    list(x = c(x[keep], end), y = c(y[keep], y[findInterval(end, x)]))
  }
  list(aj = cut_at(c(0, curves$aj$time), 100 * c(0, curves$aj$pstate[, curves$j])),
       km = cut_at(c(0, curves$km$time), 100 * c(0, 1 - curves$km$surv)))
}

# Cambia LC_CTYPE a UTF-8 si hace falta y devuelve la función que lo restaura.
use_utf8_ctype <- function() {
  if (isTRUE(l10n_info()[["UTF-8"]])) return(function() invisible(NULL))
  old <- Sys.getlocale("LC_CTYPE")
  for (loc in c("C.UTF-8", "en_US.UTF-8", "es_ES.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) {
      return(function() invisible(Sys.setlocale("LC_CTYPE", old)))
    }
  }
  function() invisible(NULL)
}


# -----------------------------------------------------------------------------
# Modelos
# -----------------------------------------------------------------------------

#' Ajusta con las mismas covariables:
#'   - Fine-Gray (sHR) del evento de interés: survival::finegray() crea el
#'     conjunto ampliado con pesos de censura y coxph() ponderado da el sHR;
#'     el error estándar es robusto, agrupado por paciente;
#'   - Cox por causa (csHR) del evento de interés y del competidor.
#' Usa casos completos en las covariables (avisa de cuántas filas quita).
#' Devuelve una fila por covariable y modelo con HR, IC 95 % y p; los ajustes
#' quedan en attr(resultado, "fits") para revisar supuestos (p. ej., cox.zph).
fit_competing_models <- function(prepared, covariates) {
  check_prepared(prepared)
  if (length(covariates) == 0) stop("Indica al menos una covariable.", call. = FALSE)
  absent <- setdiff(covariates, names(prepared))
  if (length(absent) > 0) {
    stop("Estas covariables no están en los datos preparados: ", paste(absent, collapse = ", "),
         ". Pásalas en 'covariates' de prepare_competing_data().", call. = FALSE)
  }
  d <- prepared[, c("time", "status", covariates), drop = FALSE]
  complete <- stats::complete.cases(d)
  if (!all(complete)) {
    message("Modelos: se excluyen ", sum(!complete), " filas con NA en las covariables ",
            "(casos completos, n = ", sum(complete), ").")
  }
  d <- d[complete, , drop = FALSE]
  if (!any(d$status == 1)) stop("Sin NA no queda ningún evento de interés.", call. = FALSE)
  # Una covariable con un solo valor no se puede estimar: coxph daría NA o un error
  # de contrastes. Los niveles sin pacientes de un factor se quitan.
  for (v in covariates) {
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
    if (length(unique(d[[v]])) < 2) {
      stop("La covariable '", v, "' no varía (un solo valor: ", as.character(d[[v]][1]),
           ") en los ", nrow(d), " casos completos. Quítala de 'covariates'.", call. = FALSE)
    }
  }
  d$subject_id <- seq_len(nrow(d))
  rhs <- paste(covariates, collapse = " + ")

  fg_data <- finegray(Surv(time, factor(status, 0:2, c("censored", "event", "competing"))) ~ .,
                      data = d, etype = "event")
  fits <- list(
    fine_gray = coxph(stats::as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs)),
                      data = fg_data, weights = fgwt, cluster = subject_id),
    cause_event = coxph(stats::as.formula(paste("Surv(time, status == 1) ~", rhs)), data = d)
  )
  if (any(d$status == 2)) {
    fits$cause_competing <- coxph(stats::as.formula(paste("Surv(time, status == 2) ~", rhs)),
                                  data = d)
  }

  event_label <- label_of(prepared, "event")
  competing_label <- label_of(prepared, "competing")
  rows <- list(
    tidy_cox(fits$fine_gray, "Fine-Gray", "sHR", event_label, sum(d$status == 1)),
    tidy_cox(fits$cause_event, "Cox por causa", "csHR", event_label, sum(d$status == 1))
  )
  if (!is.null(fits$cause_competing)) {
    rows[[3]] <- tidy_cox(fits$cause_competing, "Cox por causa", "csHR", competing_label,
                          sum(d$status == 2))
  }
  out <- do.call(rbind, rows)
  out$n <- nrow(d)
  no_estimate <- unique(out$term[is.na(out$hr)])
  if (length(no_estimate) > 0) {
    message("Modelos: sin estimación (NA) para ", paste(no_estimate, collapse = ", "),
            ": es combinación lineal de otras covariables.")
  }
  n_coef <- length(stats::coef(fits$cause_event))
  if (sum(d$status == 1) < 10 * n_coef) {
    message("Modelos: solo ", sum(d$status == 1), " eventos de interés para ", n_coef,
            " coeficientes (menos de 10 por coeficiente): HR inestables.")
  }
  # Orden: por covariable y, dentro de cada una, sHR, csHR del evento y csHR del competidor.
  model_index <- rep(seq_along(rows), times = vapply(rows, nrow, integer(1)))
  out <- out[order(match(out$term, unique(out$term)), model_index), ]
  rownames(out) <- NULL
  attr(out, "fits") <- fits
  out
}

tidy_cox <- function(fit, model, measure, outcome, n_events) {
  s <- summary(fit)
  data.frame(model = model, measure = measure, outcome = outcome,
             term = rownames(s$coefficients),
             hr = s$conf.int[, "exp(coef)"],
             lower = s$conf.int[, "lower .95"],
             upper = s$conf.int[, "upper .95"],
             p_value = s$coefficients[, "Pr(>|z|)"],
             n_events = n_events,
             row.names = NULL, stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Resúmenes de texto y exportación
# -----------------------------------------------------------------------------

format_num <- function(x, digits = 1) {
  ifelse(is.na(x), "-", formatC(x, format = "f", digits = digits, decimal.mark = ","))
}

# Tiempos sin notación científica (100000, no 1e+05) y con coma decimal.
format_time <- function(x) {
  vapply(x, function(v) format(v, scientific = FALSE, decimal.mark = ",", trim = TRUE),
         character(1))
}

# HR con 3 decimales; los enormes o minúsculos (modelos que no convergen), en
# notación científica para que no descuadren la tabla.
format_hr <- function(x) {
  extreme <- !is.na(x) & is.finite(x) & x != 0 & (x >= 1000 | x < 0.001)
  ifelse(extreme, formatC(x, format = "e", digits = 2, decimal.mark = ","), format_num(x, 3))
}

format_p <- function(p) {
  ifelse(is.na(p), "-", ifelse(p < 0.001, "< 0,001", format_num(p, 3)))
}

#' Líneas de texto con la tabla de compare_aj_km() en formato español.
#' Las cabeceras son ASCII para que las columnas queden alineadas en cualquier locale.
format_comparison <- function(tab) {
  ci <- paste0(format_num(100 * tab$aj, 1), " (", format_num(100 * tab$aj_lower, 1),
               " a ", format_num(100 * tab$aj_upper, 1), ")")
  no_ci <- is.na(tab$aj_lower) | is.na(tab$aj_upper)
  ci[no_ci] <- format_num(100 * tab$aj[no_ci], 1)  # "-" si tampoco hay estimación
  note <- ifelse(tab$n_risk == 0, "  sin seguimiento",
                 ifelse(tab$few_at_risk, "  pocos en riesgo", ""))
  diff <- ifelse(is.na(tab$diff_pp), "-", paste0("+", format_num(tab$diff_pp, 1)))
  c(sprintf("%8s %9s %20s %8s %8s %8s", "tiempo", "en riesgo", "AJ % (IC 95 %)",
            "1-KM %", "dif. pp", "1-KM/AJ"),
    paste0(sprintf("%8s %9d %20s %8s %8s %8s", format_time(tab$time),
                   as.integer(tab$n_risk), ci, format_num(100 * tab$one_minus_km, 1),
                   diff, format_num(tab$ratio, 2)), note))
}

#' Líneas de texto con la tabla de fit_competing_models() en formato español.
format_models <- function(models) {
  hr <- paste0(format_hr(models$hr), " (", format_hr(models$lower), " a ",
               format_hr(models$upper), ")")
  c(sprintf("%-10s %-5s %-14s %-24s %-8s %s", "variable", "HR", "modelo",
            "HR (IC 95 %)", "p", "resultado"),
    sprintf("%-10s %-5s %-14s %-24s %-8s %s", models$term, models$measure, models$model,
            hr, format_p(models$p_value), models$outcome))
}

#' Escribe una tabla CSV con 6 cifras significativas. `sep` y `dec` permiten el
#' formato de Excel en español (sep = ";", dec = ",").
write_table <- function(df, file, sep = ",", dec = ".") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = 6)
  utils::write.table(df, file, sep = sep, dec = dec, row.names = FALSE, qmethod = "double")
  invisible(file)
}

#' Escribe los metadatos de un análisis (formato "clave: valor", legible con
#' read.dcf()): versión, fecha, R, survival, entrada, recuentos y parámetros.
write_metadata <- function(file, prepared, input, parameters = list()) {
  check_prepared(prepared)
  fields <- c(
    herramienta = paste("riesgos-competitivos", rc_version),
    fecha = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    version_R = R.version.string,
    version_survival = utils::packageDescription("survival")$Version,
    entrada = input,
    n = nrow(prepared),
    n_evento = sum(prepared$status == 1),
    n_competidor = sum(prepared$status == 2),
    n_censura = sum(prepared$status == 0),
    etiqueta_evento = label_of(prepared, "event"),
    etiqueta_competidor = label_of(prepared, "competing"),
    vapply(parameters, function(x) paste(x, collapse = ", "), character(1)),
    aviso = "Uso docente y de investigación; no apto para decisiones clínicas."
  )
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  # useBytes evita que las tildes salgan escapadas en sesiones sin UTF-8.
  writeLines(paste0(names(fields), ": ", gsub("\n", " ", fields)), file, useBytes = TRUE)
  invisible(file)
}
