# =============================================================================
# Tests sin frameworks: solo base R + survival (cmprsk es opcional).
# Ejecuta desde la carpeta del repositorio:   Rscript tests/run_tests.R
# Sale con código 1 si falla alguna comprobación. Con la variable de entorno
# REQUIRE_CMPRSK=true (así en CI), la falta de cmprsk cuenta como fallo.
# =============================================================================

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
# En Linux y Mac, Rscript codifica como "~+~" los espacios de la ruta.
script <- gsub("~+~", " ", script, fixed = TRUE)
if (length(script) == 1) setwd(file.path(dirname(normalizePath(script)), ".."))
if (!file.exists("R/riesgos_competitivos.R")) {
  stop("No encuentro R/riesgos_competitivos.R. Ejecuta los tests desde la carpeta del repositorio.")
}
source("R/riesgos_competitivos.R")

# ---- Mini arnés de pruebas --------------------------------------------------
passed <- 0L
failed <- character(0)
skipped <- character(0)

check <- function(description, condition) {
  if (isTRUE(condition)) {
    passed <<- passed + 1L
    cat("  ok       ", description, "\n")
  } else {
    failed <<- c(failed, description)
    cat("  FALLO    ", description, "\n")
  }
}

skip <- function(description, reason) {
  skipped <<- c(skipped, description)
  cat("  OMITIDA  ", description, "--", reason, "\n")
}

# Ejecuta un bloque; un error inesperado cuenta como fallo y no para el resto.
section <- function(title, code) {
  cat("\n", title, "\n", sep = "")
  tryCatch(code, error = function(e) {
    failed <<- c(failed, paste0(title, ": error inesperado"))
    cat("  FALLO     error inesperado:", conditionMessage(e), "\n")
  })
}

error_message <- function(expr) {
  tryCatch({
    suppressMessages(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

check_error <- function(description, expr, pattern) {
  msg <- error_message(expr)
  check(paste0(description, if (is.na(msg)) " (no dio error)" else ""),
        !is.na(msg) && grepl(pattern, msg))
}

quiet <- function(expr) suppressMessages(expr)

# ---- Datos --------------------------------------------------------------------
# mgus2 como en la viñeta "compete" de survival (tiempo en meses).
mgus <- survival::mgus2
mgus$etime <- ifelse(mgus$pstat == 0, mgus$futime, mgus$ptime)
mgus$event <- ifelse(mgus$pstat == 0, 2 * mgus$death, 1)
prep_mgus <- function(event_codes = 1, competing_codes = 2, data = mgus) {
  quiet(prepare_competing_data(data, "etime", "event", event_codes = event_codes,
                               competing_codes = competing_codes, censor_codes = 0,
                               covariates = c("age", "sex"),
                               event_label = "PCM", competing_label = "muerte"))
}
d_mgus <- prep_mgus()

# Sintéticos con semilla fija; los tiempos redondeados crean muchos empates.
set.seed(20260930)
sim <- data.frame(time = round(stats::rexp(400, 1 / 8), 1),
                  status = sample(0:2, 400, replace = TRUE, prob = c(0.3, 0.3, 0.4)))
prep_sim <- function(data = sim, event_codes = 1, competing_codes = 2) {
  quiet(prepare_competing_data(data, "time", "status", event_codes = event_codes,
                               competing_codes = competing_codes, censor_codes = 0))
}
d_sim <- prep_sim()

# Rejilla: todos los tiempos observados y puntos intermedios, dentro del seguimiento.
time_grid <- function(d) {
  g <- sort(unique(c(0, d$time, d$time + 0.5)))
  g[g <= max(d$time)]
}

# Aalen-Johansen escrito a mano, independiente de survival:
#   F_k(t) = suma, en los tiempos de evento t_j <= t, de S(t_j-) * d_kj / n_j,
#   con S el Kaplan-Meier de "cualquier evento".
aj_by_hand <- function(time, status, cause) {
  event_times <- sort(unique(time[status != 0]))
  surv <- 1
  cif <- 0
  values <- numeric(length(event_times))
  for (i in seq_along(event_times)) {
    at_t <- time == event_times[i]
    n_risk <- sum(time >= event_times[i])
    cif <- cif + surv * sum(at_t & status == cause) / n_risk
    surv <- surv * (1 - sum(at_t & status != 0) / n_risk)
    values[i] <- cif
  }
  function(t) c(0, values)[findInterval(t, event_times) + 1]
}

# ---- 1. Aalen-Johansen frente a la implementación a mano --------------------
section("1. AJ de survfit (vía compare_aj_km) igual que el escrito a mano", {
  for (case in list(list("mgus2", d_mgus), list("sintéticos con empates", d_sim))) {
    d <- case[[2]]
    grid <- time_grid(d)
    tab <- quiet(compare_aj_km(d, grid))
    diff <- max(abs(tab$aj - aj_by_hand(d$time, d$status, 1)(grid)))
    check(sprintf("%s: %d tiempos, diferencia máxima %.1e (tolerancia 1e-10)",
                  case[[1]], length(grid), diff), diff < 1e-10)
  }
})

# ---- 2. Sin competidores, AJ == 1 - KM --------------------------------------
section("2. Sin eventos competidores, AJ == 1 - KM", {
  no_comp <- sim
  no_comp$status[no_comp$status == 2] <- 0
  msgs <- character(0)
  d0 <- withCallingHandlers(
    prepare_competing_data(no_comp, "time", "status", 1, 2, 0),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  check("avisa de que no hay competidores", any(grepl("No hay eventos competidores", msgs)))
  tab <- quiet(compare_aj_km(d0, time_grid(d0)))
  diff <- max(abs(tab$aj - tab$one_minus_km))
  check(sprintf("diferencia máxima %.1e (tolerancia 1e-12)", diff), diff < 1e-12)
})

# ---- 3. Con competidores, 1 - KM >= AJ --------------------------------------
section("3. Con competidores, 1 - KM >= AJ en todos los tiempos", {
  for (case in list(list("mgus2", d_mgus), list("sintéticos", d_sim))) {
    tab <- quiet(compare_aj_km(case[[2]], time_grid(case[[2]])))
    check(sprintf("%s: 1 - KM >= AJ en %d tiempos y > al final", case[[1]], nrow(tab)),
          all(tab$one_minus_km - tab$aj >= -1e-12) &&
            tab$one_minus_km[nrow(tab)] > tab$aj[nrow(tab)])
  }
})

# ---- 4. Incidencias de todas las causas + supervivencia global = 1 ----------
section("4. Incidencia de cada causa + supervivencia global = 1", {
  cases <- list(list("mgus2", prep_mgus(1, 2), prep_mgus(2, 1)),
                list("sintéticos", prep_sim(sim, 1, 2), prep_sim(sim, 2, 1)))
  for (case in cases) {
    grid <- time_grid(case[[2]])
    f1 <- quiet(compare_aj_km(case[[2]], grid))$aj
    f2 <- quiet(compare_aj_km(case[[3]], grid))$aj
    overall <- summary(survfit(Surv(time, status != 0) ~ 1, data = case[[2]]),
                       times = grid, extend = TRUE)$surv
    diff <- max(abs(f1 + f2 + overall - 1))
    check(sprintf("%s: diferencia máxima %.1e (tolerancia 1e-10)", case[[1]], diff), diff < 1e-10)
  }
})

# ---- 5. Cifras de mgus2 frente a la referencia versionada -------------------
section("5. mgus2 frente a tests/referencia_mgus2*.csv (survival 3.5-8)", {
  ref <- utils::read.csv("tests/referencia_mgus2.csv", comment.char = "#")
  tab <- quiet(compare_aj_km(d_mgus, ref$time))
  check("pacientes en riesgo idénticos", identical(as.numeric(tab$n_risk), as.numeric(ref$n_risk)))
  check("AJ y 1 - KM (tolerancia 1e-8)",
        max(abs(tab$aj - ref$aj), abs(tab$one_minus_km - ref$one_minus_km)) < 1e-8)
  # Los IC dependen del método de varianza de survival: tolerancia de 0,01 puntos.
  check("IC 95 % de AJ (tolerancia 1e-4)",
        max(abs(tab$aj_lower - ref$aj_lower), abs(tab$aj_upper - ref$aj_upper)) < 1e-4)
  check("pocos en riesgo (< 20) solo a 300 y 360 meses",
        identical(tab$time[tab$few_at_risk], c(300, 360)))
  # Cifras de control calculadas aparte, redondeadas a 4 decimales.
  ctl <- tab[tab$time %in% c(180, 240), ]
  check("control: 15 años AJ 0,0871, 1 - KM 0,1596, 178 en riesgo; 20 años 0,0998, 0,2096, 57",
        identical(round(ctl$aj, 4), c(0.0871, 0.0998)) &&
          identical(round(ctl$one_minus_km, 4), c(0.1596, 0.2096)) &&
          identical(as.numeric(ctl$n_risk), c(178, 57)))
  check("diferencia y razón: +7,26 y +10,97 puntos; x1,83 y x2,10",
        identical(round(ctl$diff_pp, 2), c(7.26, 10.97)) &&
          identical(round(ctl$ratio, 2), c(1.83, 2.10)))

  refm <- utils::read.csv("tests/referencia_mgus2_modelos.csv", comment.char = "#",
                          stringsAsFactors = FALSE)
  mod <- quiet(fit_competing_models(d_mgus, c("age", "sex")))
  m <- match(paste(refm$model, refm$outcome, refm$term), paste(mod$model, mod$outcome, mod$term))
  check("modelos: mismas filas que la referencia", !anyNA(m) && nrow(mod) == nrow(refm))
  check("modelos: HR (tolerancia relativa 1e-6)", max(abs(mod$hr[m] / refm$hr - 1)) < 1e-6)
  check("modelos: IC 95 % (tolerancia relativa 1e-4)",
        max(abs(mod$lower[m] / refm$lower - 1), abs(mod$upper[m] / refm$upper - 1)) < 1e-4)
  check("modelos: p (tolerancia relativa 1e-4)",
        isTRUE(all.equal(mod$p_value[m], refm$p_value, tolerance = 1e-4)))
  fg <- mod[mod$model == "Fine-Gray", ]
  check("control: Fine-Gray log sHR edad -0,01730 y sexo M -0,25976",
        identical(round(log(fg$hr), 5), c(-0.01730, -0.25976)))
  check("modelos ordenados por covariable: sHR, csHR evento, csHR competidor",
        identical(mod$term, rep(c("age", "sexM"), each = 3)) &&
          identical(mod$measure, rep(c("sHR", "csHR", "csHR"), 2)) &&
          all(mod$lower < mod$hr & mod$hr < mod$upper))
})

# ---- 6. Fine-Gray de survival frente a cmprsk::crr --------------------------
section("6. Fine-Gray (finegray + coxph) frente a cmprsk::crr", {
  if (requireNamespace("cmprsk", quietly = TRUE)) {
    # a) Sin empates, las dos implementaciones dan los mismos coeficientes. crr para
    #    por defecto con gtol = 1e-6 (diferencias ~1e-6); con 1e-9 converge del todo.
    set.seed(7)
    n <- 500
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1, 0.5)
    t1 <- stats::rexp(n, 0.10 * exp(0.5 * x1))
    t2 <- stats::rexp(n, 0.15 * exp(0.3 * x1 - 0.4 * x2))
    cens <- stats::rexp(n, 0.05)
    tt <- pmin(t1, t2, cens)
    st <- ifelse(tt == t1, 1, ifelse(tt == t2, 2, 0))
    d <- quiet(prepare_competing_data(data.frame(tt, st, x1, x2), "tt", "st", 1, 2, 0,
                                      covariates = c("x1", "x2")))
    mod <- quiet(fit_competing_models(d, c("x1", "x2")))
    crr <- cmprsk::crr(tt, st, cov1 = cbind(x1, x2), failcode = 1, cencode = 0, gtol = 1e-9)
    diff <- max(abs(log(mod$hr[mod$model == "Fine-Gray"]) - crr$coef))
    check(sprintf("sintéticos sin empates: diferencia máxima %.1e (tolerancia 1e-8)", diff),
          diff < 1e-8)

    # b) mgus2 tiene empates: difieren en el 4.º decimal, muy por debajo del error estándar.
    crr <- cmprsk::crr(d_mgus$time, d_mgus$status, failcode = 1, cencode = 0, gtol = 1e-9,
                       cov1 = cbind(age = d_mgus$age, sexM = as.numeric(d_mgus$sex == "M")))
    fg <- quiet(fit_competing_models(d_mgus, c("age", "sex")))
    fg <- fg[fg$model == "Fine-Gray", ]
    se_crr <- sqrt(diag(crr$var))
    se_fg <- (log(fg$upper) - log(fg$lower)) / (2 * stats::qnorm(0.975))
    check("mgus2: |diferencia de coeficientes| < 5 % del error estándar de crr",
          all(abs(log(fg$hr) - crr$coef) < 0.05 * se_crr))
    check("mgus2: error estándar robusto a menos del 10 % del de crr",
          all(abs(se_fg / se_crr - 1) < 0.10))
  } else if (identical(Sys.getenv("REQUIRE_CMPRSK"), "true")) {
    check("cmprsk instalado (REQUIRE_CMPRSK=true)", FALSE)
  } else {
    skip("Fine-Gray frente a cmprsk::crr", "falta el paquete cmprsk: install.packages(\"cmprsk\")")
  }
})

# ---- 7. Validación de la entrada y tiempos fuera del seguimiento ------------
section("7. Validación de la entrada", {
  base <- data.frame(t = c(5, 3, 8, 2, 7), s = c(1, 0, 2, 1, 0))
  prep <- function(df, ev = 1, comp = 2, cens = 0) {
    prepare_competing_data(df, "t", "s", event_codes = ev, competing_codes = comp,
                           censor_codes = cens)
  }
  bad_t <- function(v) transform(base, t = v)
  bad_s <- function(v) transform(base, s = v)
  check_error("NA en el tiempo", prep(bad_t(c(5, NA, 8, 2, 7))), "Hay NA en 't'")
  check_error("NA en el estado", prep(bad_s(c(1, NA, 2, 1, 0))), "Hay NA en 's'")
  check_error("tiempos negativos", prep(bad_t(c(5, -1, 8, 2, 7))), "negativos en 't'")
  check_error("tiempos infinitos", prep(bad_t(c(5, Inf, 8, 2, 7))), "infinitos")
  check_error("tiempo no numérico", prep(bad_t(c("5", "3", "8", "2", "7"))), "no es num")
  check_error("códigos no declarados", prep(bad_s(c(1, 0, 2, 3, 0))), "no declarados.*'3'")
  check_error("solape entre códigos", prep(base, comp = c(2, 0)), "repetidos.*: 0")
  check_error("cero eventos de interés", prep(bad_s(c(2, 0, 2, 0, 0))), "evento de inter")
  check_error("lista de códigos vacía", prep(base, cens = character(0)), "al menos un c")
  check_error("columna inexistente", prepare_competing_data(base, "tiempo", "s", 1, 2, 0),
              "No existen.*tiempo")
  check_error("covariable con nombre reservado",
              prepare_competing_data(transform(base, status = 1), "t", "s", 1, 2, 0,
                                     covariates = "status"), "Renombra")
  check_error("datos sin preparar", compare_aj_km(base, 5), "prepare_competing_data")
  check_error("tiempos de la tabla negativos", compare_aj_km(prep(base), -1), "'times'")
  # Con Inf, survival 3.5-8 devolvía sin error la incidencia de otro estado.
  check_error("tiempos de la tabla infinitos", compare_aj_km(d_mgus, c(120, Inf)), "finitos")
  for (v in c("fgstart", "fgstop", "fgstatus", "fgwt", "subject_id", "time")) {
    df <- base
    df[[v]] <- 1
    check_error(paste0("covariable con nombre reservado: ", v),
                prepare_competing_data(df, "t", "s", 1, 2, 0, covariates = v), "Renombra")
  }

  # Códigos de texto (con espacios sobrantes) dan lo mismo que los numéricos.
  txt <- mgus
  txt$event <- c(" censura", "PCM ", "muerte")[mgus$event + 1]
  d_txt <- quiet(prepare_competing_data(txt, "etime", "event", "PCM", "muerte", "censura"))
  grid <- c(12, 60, 180, 240)
  check("códigos de texto == códigos numéricos",
        isTRUE(all.equal(quiet(compare_aj_km(d_txt, grid)), quiet(compare_aj_km(d_mgus, grid)))))

  check("tiempos desordenados o repetidos: se ordenan sin mezclar valores",
        isTRUE(all.equal(quiet(compare_aj_km(d_mgus, c(240, 180, 240))),
                         quiet(compare_aj_km(d_mgus, c(180, 240))))))

  last <- max(d_mgus$time)
  msgs <- character(0)
  tab <- withCallingHandlers(
    compare_aj_km(d_mgus, c(last, last + 1, 600)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  check("en el último tiempo de seguimiento sí hay estimación",
        !is.na(tab$aj[1]) && !is.na(tab$one_minus_km[1]) && tab$n_risk[1] > 0)
  check("más allá del seguimiento: NA, sin arrastrar el último valor",
        all(is.na(unlist(tab[2:3, c("aj", "aj_lower", "aj_upper", "one_minus_km",
                                    "diff_pp", "ratio")]))) && all(tab$n_risk[2:3] == 0))
  check("avisa de pocos en riesgo y de tiempos sin seguimiento",
        any(grepl("Pocos pacientes en riesgo", msgs)) && any(grepl("Nadie sigue", msgs)))

  d_na <- d_mgus
  d_na$age[1:10] <- NA
  mod_na <- quiet(fit_competing_models(d_na, c("age", "sex")))
  check("modelos con NA en covariables: casos completos (n = 1374)", all(mod_na$n == 1374))
})

# ---- 8. Gráfico ------------------------------------------------------------------
section("8. Gráfico", {
  f <- tempfile(fileext = ".png")
  devices <- length(grDevices::dev.list())
  ctype <- Sys.getlocale("LC_CTYPE")
  plot_aj_vs_km(d_mgus, file = f, highlight_time = 240, caption = "Prueba con tildes: años")
  png_signature <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  check("escribe un PNG no vacío",
        file.exists(f) && file.size(f) > 10000 && identical(readBin(f, "raw", 8), png_signature))
  check("cierra su dispositivo y deja el locale como estaba",
        length(grDevices::dev.list()) == devices && identical(Sys.getlocale("LC_CTYPE"), ctype))
  unlink(f)

  # Las curvas que se dibujan: cortadas en xmax y, como la tabla, en el último seguimiento.
  at240 <- quiet(compare_aj_km(d_mgus, 240))
  xy <- curve_points(d_mgus, 240)
  n_aj <- length(xy$aj$x)
  n_km <- length(xy$km$x)
  check("curvas cortadas en xmax con el valor de la tabla",
        max(xy$aj$x) == 240 && max(xy$km$x) == 240 && !is.unsorted(xy$aj$x) &&
          abs(xy$aj$y[n_aj] - 100 * at240$aj) < 1e-10 &&
          abs(xy$km$y[n_km] - 100 * at240$one_minus_km) < 1e-10)
  last <- max(d_mgus$time)
  xy <- curve_points(d_mgus, last + 120)
  check("con xmax mayor que el seguimiento, las curvas acaban en el último seguimiento",
        max(xy$aj$x) == last && max(xy$km$x) == last)
  check("curvas desde 0 y no decrecientes",
        xy$aj$x[1] == 0 && xy$aj$y[1] == 0 && xy$km$y[1] == 0 &&
          !is.unsorted(xy$aj$x) && !is.unsorted(xy$aj$y) && !is.unsorted(xy$km$y))

  # Un PNG con ejes y leyenda pero sin curvas también es "un PNG no vacío": en SVG
  # se puede comprobar que las dos curvas están dibujadas, cada una con su color.
  # svg() necesita cairo; sin él (p. ej., macOS sin XQuartz) solo avisa y no abre nada.
  f <- tempfile(fileext = ".svg")
  before <- grDevices::dev.cur()
  suppressWarnings(grDevices::svg(f, width = 8, height = 5.4))
  if (!identical(grDevices::dev.cur(), before)) {
    plot_aj_vs_km(d_mgus, highlight_time = 240)
    grDevices::dev.off()
    svg <- readLines(f, warn = FALSE)
    # Azul #0072B2 (AJ) y naranja #D55E00 (1 - KM), en % como los escribe cairo. La
    # leyenda y los puntos dan trazos cortos; la curva es el más largo de cada color.
    curve_of <- function(colour) {
      hits <- grep(colour, svg, value = TRUE)
      if (length(hits) == 0) return(list(chars = 0L, top = NA_real_))
      path <- sub('.* d="([^"]*)".*', "\\1", hits[which.max(nchar(hits))])
      xy <- as.numeric(regmatches(path, gregexpr("[0-9.]+", path))[[1]])
      # En SVG la y crece hacia abajo: el punto más alto de la curva es la y mínima.
      list(chars = nchar(path), top = min(xy[c(FALSE, TRUE)]))
    }
    aj_curve <- curve_of("44[.]70")
    km_curve <- curve_of("83[.]52")
    check(sprintf("SVG: curva AJ dibujada (%d caracteres) y curva 1 - KM (%d)",
                  aj_curve$chars, km_curve$chars),
          aj_curve$chars > 1500 && km_curve$chars > 1500)
    check("SVG: la curva naranja (1 - KM) acaba por encima de la azul (AJ)",
          isTRUE(km_curve$top < aj_curve$top - 10))
    unlink(f)
  } else {
    skip("curvas dibujadas (SVG)", "esta instalación de R no puede abrir svg() (falta cairo)")
  }
})

# ---- 9. Recuentos, umbral de pocos en riesgo y avisos de los modelos ---------
section("9. Recuentos, umbral de pocos en riesgo y avisos de los modelos", {
  mod <- quiet(fit_competing_models(d_mgus, c("age", "sex")))
  check("n_events: 115 PCM y 860 muertes; n = 1384",
        all(mod$n_events[mod$outcome == "PCM"] == sum(d_mgus$status == 1)) &&
          all(mod$n_events[mod$outcome == "muerte"] == sum(d_mgus$status == 2)) &&
          sum(d_mgus$status == 1) == 115 && sum(d_mgus$status == 2) == 860 &&
          all(mod$n == 1384))
  # A 240 meses quedan 57 en riesgo: con el umbral en 57 no es "pocos"; con 58, sí.
  check("pocos en riesgo es estrictamente menor que min_at_risk",
        !quiet(compare_aj_km(d_mgus, 240, min_at_risk = 57))$few_at_risk &&
          quiet(compare_aj_km(d_mgus, 240, min_at_risk = 58))$few_at_risk)

  messages_of <- function(expr) {
    msgs <- character(0)
    suppressWarnings(withCallingHandlers(expr, message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
    msgs
  }
  d_alias <- d_mgus
  d_alias$age2 <- 2 * d_alias$age
  check("avisa de las covariables sin estimación (combinación lineal de otras)",
        any(grepl("sin estimación \\(NA\\) para age2",
                  messages_of(fit_competing_models(d_alias, c("age", "age2"))))))
  few <- d_mgus[c(which(d_mgus$status == 1)[1:6], which(d_mgus$status == 2)[1:20],
                  which(d_mgus$status == 0)[1:20]), ]
  check("avisa de pocos eventos por coeficiente",
        any(grepl("solo 6 eventos de interés para 2 coeficientes",
                  messages_of(fit_competing_models(few, c("age", "sex"))))))
  check("con eventos suficientes no avisa",
        length(messages_of(fit_competing_models(d_mgus, c("age", "sex")))) == 0)
})

# ---- 9b. Casos límite de los modelos y de los parámetros ----------------------
section("9b. Casos límite de los modelos y de los parámetros", {
  d_const <- d_mgus
  d_const$k <- 1
  check_error("covariable numérica constante", fit_competing_models(d_const, c("age", "k")),
              "'k' no varía")
  d_const$k <- "F"
  check_error("covariable de texto con un solo valor",
              fit_competing_models(d_const, c("age", "k")), "'k' no varía")
  check_error("factor con un solo grupo en el subconjunto",
              fit_competing_models(d_mgus[d_mgus$sex == "F", ], c("age", "sex")), "'sex' no varía")
  d_lev <- d_mgus
  d_lev$sex <- factor(d_lev$sex, levels = c("F", "M", "otro"))
  mod <- quiet(fit_competing_models(d_lev, c("age", "sex")))
  check("factor con un nivel sin pacientes: se quita el nivel, sin filas NA",
        !anyNA(mod$hr) && identical(unique(mod$term), c("age", "sexM")))

  inf <- transform(mgus, x = ifelse(seq_len(nrow(mgus)) == 3, Inf, age))
  check_error("Inf en una covariable",
              prepare_competing_data(inf, "etime", "event", 1, 2, 0, covariates = "x"),
              "infinitos en 'x' \\(filas: 3\\)")

  # Sin competidores: dos modelos, y Fine-Gray da el mismo HR que Cox.
  no_comp <- mgus
  no_comp$event[no_comp$event == 2] <- 0
  d0 <- quiet(prepare_competing_data(no_comp, "etime", "event", 1, 2, 0,
                                     covariates = c("age", "sex")))
  mod <- quiet(fit_competing_models(d0, c("age", "sex")))
  check("sin competidores: solo Fine-Gray y Cox del evento, con el mismo HR",
        nrow(mod) == 4 && all(mod$outcome == "evento") &&
          max(abs(mod$hr[mod$model == "Fine-Gray"] / mod$hr[mod$model == "Cox por causa"] - 1)) < 1e-8)

  for (v in list(NA, -5, "a", c(10, 20), Inf)) {
    check_error(paste0("min_at_risk no válido: ", paste(deparse(v), collapse = "")),
                compare_aj_km(d_mgus, 240, min_at_risk = v), "'min_at_risk'")
  }
  check_error("min_at_risk no válido en el gráfico",
              plot_aj_vs_km(d_mgus, file = tempfile(fileext = ".png"), min_at_risk = NA),
              "'min_at_risk'")
  check_error("xmax no válido en el gráfico",
              plot_aj_vs_km(d_mgus, file = tempfile(fileext = ".png"), xmax = Inf), "'xmax'")
})

# ---- 10. Formato y exportación ------------------------------------------------
section("10. Formato y exportación", {
  check("format_num: coma decimal y '-' para NA",
        identical(format_num(c(1234.5, 0.25, NA), 2), c("1234,50", "0,25", "-")))
  check("format_p: '< 0,001' solo por debajo de 0,001",
        identical(format_p(c(0.0004, 0.001, 0.0124, 0.5, NA)),
                  c("< 0,001", "0,001", "0,012", "0,500", "-")))
  check("format_time: sin notación científica",
        identical(format_time(c(0.5, 60, 100000)), c("0,5", "60", "100000")))
  check("format_hr: notación científica para HR enormes",
        identical(format_hr(c(1.5, 1e21, NA)), c("1,500", "1,00e+21", "-")))

  tab <- quiet(compare_aj_km(d_mgus, c(180, 240)))
  lines <- format_comparison(tab)
  check("format_comparison: fila de 240 meses",
        length(lines) == 3 &&
          grepl("^ +240 +57 +10,0 \\(8,2 a 12,1\\) +21,0 +\\+11,0 +2,10$", lines[3]))
  mod <- quiet(fit_competing_models(d_mgus, c("age", "sex")))
  lines <- format_models(mod)
  check("format_models: fila de la edad en Fine-Gray",
        length(lines) == 7 &&
          grepl("^age +sHR +Fine-Gray +0,983 \\(0,972 a 0,994\\) +0,002 +PCM$", lines[2]))

  # CSV: 6 cifras significativas, y el formato de Excel en español se relee igual.
  df <- data.frame(a = c(0.123456789, 123456.789), b = c("x", "y"), n = 1:2)
  f <- tempfile(fileext = ".csv")
  write_table(df, f)
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  check("write_table: 6 cifras significativas",
        identical(back$a, c(0.123457, 123457)) && identical(back$b, c("x", "y")) &&
          identical(back$n, 1:2))
  write_table(tab, f, sep = ";", dec = ",")
  back <- utils::read.csv2(f)
  check("write_table con ';' y coma decimal: la tabla se relee igual (6 cifras)",
        identical(names(back), names(tab)) &&
          isTRUE(all.equal(back$aj, signif(tab$aj, 6), tolerance = 1e-12)) &&
          isTRUE(all.equal(back$ratio, signif(tab$ratio, 6), tolerance = 1e-12)) &&
          identical(back$n_risk, c(178L, 57L)))
  unlink(f)

  f <- tempfile(fileext = ".txt")
  write_metadata(f, d_mgus, input = "prueba", parameters = list(tiempos = c(60, 120)))
  meta <- read.dcf(f)[1, ]
  check("write_metadata: recuentos, entrada, parámetros y versiones",
        identical(unname(meta[c("n", "n_evento", "n_competidor", "n_censura")]),
                  c("1384", "115", "860", "409")) &&
          meta[["entrada"]] == "prueba" && meta[["tiempos"]] == "60, 120" &&
          meta[["etiqueta_evento"]] == "PCM" &&
          meta[["version_survival"]] == utils::packageDescription("survival")$Version &&
          grepl("^riesgos-competitivos ", meta[["herramienta"]]))
  unlink(f)
})

# ---- 11. Los scripts, de principio a fin ----------------------------------------
# Se ejecutan en una copia dentro de una carpeta con espacios y llamándolos desde
# otra carpeta: es el caso que fallaba en Linux y Mac (Rscript codifica los
# espacios de la ruta como "~+~").
section("11. plantilla.R y ejemplo_mgus2.R en una carpeta con espacios", {
  copy <- file.path(tempfile("rc "), "mi carpeta")
  dir.create(file.path(copy, "R"), recursive = TRUE)
  file.copy(c("plantilla.R", "ejemplo_mgus2.R"), copy)
  file.copy("R/riesgos_competitivos.R", file.path(copy, "R"))
  rscript <- file.path(R.home("bin"), "Rscript")
  run <- function(script) {
    out <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(file.path(copy, script))),
                                    stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) cat(paste0("      ", utils::tail(out, 5), "\n"), sep = "")
    is.null(status) || status == 0
  }

  check("plantilla.R termina sin error", run("plantilla.R"))
  out <- file.path(copy, "resultados")
  tab_file <- file.path(out, "plantilla_tabla_aj_km.csv")
  expected <- quiet(compare_aj_km(d_mgus, c(60, 120, 180, 240)))
  check("plantilla.R: su tabla coincide con compare_aj_km() (6 cifras)",
        file.exists(tab_file) && {
          got <- utils::read.csv(tab_file)
          isTRUE(all.equal(got$aj, signif(expected$aj, 6), tolerance = 1e-12)) &&
            isTRUE(all.equal(got$one_minus_km, signif(expected$one_minus_km, 6), tolerance = 1e-12)) &&
            identical(as.numeric(got$n_risk), as.numeric(expected$n_risk))
        })
  meta_file <- file.path(out, "plantilla_metadatos.txt")
  check("plantilla.R: metadatos con el md5 de la entrada y los recuentos",
        file.exists(meta_file) && {
          meta <- read.dcf(meta_file)[1, ]
          grepl("md5 [0-9a-f]{32}", meta[["entrada"]]) && meta[["n_evento"]] == "115" &&
            meta[["n_competidor"]] == "860"
        })
  check("plantilla.R: escribe figura y modelos",
        all(file.exists(file.path(out, c("plantilla_figura.png", "plantilla_modelos.csv")))))

  check("ejemplo_mgus2.R termina sin error", run("ejemplo_mgus2.R"))
  ex_file <- file.path(out, "mgus2_tabla_aj_km.csv")
  check("ejemplo_mgus2.R: a 20 años, AJ 0,0998 y 1 - KM 0,2096",
        file.exists(ex_file) && {
          got <- utils::read.csv(ex_file)
          at20 <- got[got$time == 20, ]
          nrow(at20) == 1 && round(at20$aj, 4) == 0.0998 && round(at20$one_minus_km, 4) == 0.2096
        })
  unlink(dirname(copy), recursive = TRUE)
})

# ---- Resumen -------------------------------------------------------------------
total <- passed + length(failed)
cat("\nResumen: ", total, " comprobaciones, ", passed, " correctas, ", length(failed),
    " fallidas, ", length(skipped), " omitidas.\n", sep = "")
if (length(skipped) > 0) cat("AVISO, omitidas:", paste(skipped, collapse = "; "), "\n")
if (length(failed) > 0) {
  cat("Fallidas:\n", paste0("  - ", failed, "\n"), sep = "")
  if (interactive()) stop("Hay tests fallidos.") else quit(status = 1)
}
