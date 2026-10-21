# Exact output checks use asymmetric draws, without fitting a model.
interval_fixture <- function(n_groups = 2, covariates = TRUE) {
  z <- seq(-1, 2, length.out = 101)^3 / 4
  fit <- cbind(1, seq(-0.8, 0.6, length.out = 101))
  attr(fit, "draws") <- list(theta = z + 2, sigma = exp(z)/5,
                             theta_k = outer(z, seq_len(n_groups), "+"),
                             beta = matrix(z/2, ncol = 1),
                             omega = matrix(exp(z)/4, ncol = 1))
  inputs <- structure(list(X = matrix(0, n_groups, 1), M = 1, c = 1.96),
                      n_groups = n_groups, group_label = paste0("Study ", seq_len(n_groups)),
                      covariate_coding = "x")
  dat <- data.frame(group = attr(inputs, "group_label"),
                    tau = seq_len(n_groups), se = seq_len(n_groups)/2)
  structure(list(model = "rubin", fit = fit, inputs = inputs,
                 effects = "mean", n_parameters = 1, n_groups = n_groups,
                 pooling = "partial", covariates = if(covariates) "x" else character(),
                 data = dat, summary_data = dat), class = "baggr", ppd = FALSE)
}

mock_interval_extract <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    extract = function(object, pars, ...) list(attr(object, "draws")[[pars]]),
    .package = "rstan", .env = env)
}

expected_interval <- function(x, prob) {
  c(mean = mean(x), median = median(x), sd = sd(x),
    lower = unname(quantile(x, (1 - prob)/2)),
    upper = unname(quantile(x, (1 + prob)/2)), width = prob)
}

expect_printed_table <- function(output, tab) {
  expect_true(all(capture.output(print(tab, digits = 7)) %in% output))
}

test_that("prob means contained mass and legacy arguments remain compatible", {
  mock_interval_extract()
  bg <- interval_fixture()
  z <- attr(bg$fit, "draws")$theta
  for(p in c(0, .05, .5, .8, .95, 1)) {
    got <- mint(z, prob = p, median = TRUE, sd = TRUE)
    expect_equal(unname(got), unname(expected_interval(z, p)))
    expect_identical(mint(z, p), mint(z, int = p))
    for(f in list(group_effects, fixed_effects, treatment_effect,
                  hypermean, hypersd, selection, vivs_cor, mutau_cor)) {
      expect_equal(f(bg, prob = p), f(bg, interval = p))
      expect_error(f(bg, prob = p, interval = p), "Supply only one")
    }
    set.seed(130)
    a <- predict(bg, summary = TRUE, prob = p)
    set.seed(130)
    expect_equal(a, predict(bg, summary = TRUE, interval = p))
    expect_equal(group_effects(bg, TRUE, NULL, p),
                 group_effects(bg, summary = TRUE, prob = p))
    expect_equal(fixed_effects(bg, TRUE, NULL, p),
                 fixed_effects(bg, summary = TRUE, prob = p))
    expect_equal(treatment_effect(bg, TRUE, NULL, p),
                 treatment_effect(bg, summary = TRUE, prob = p))
  }
  for(p in list(-.1, 1.1, 95, NA_real_, Inf, numeric(), c(.5,.95), "0.95")) {
    expect_error(mint(z, prob = p), "single finite number")
    expect_error(group_effects(bg, prob = p), "single finite number")
  }
  expect_error(mint(z, prob = .5, int = .5), "Supply only one")
  expect_error(predict(bg, prob = .5, interval = .5), "Supply only one")
})

test_that("summary values describe the requested draws and transformation", {
  mock_interval_extract()
  bg <- interval_fixture()
  draws <- attr(bg$fit, "draws")
  for(p in c(.05, .5, .8, .95)) {
    for(tr in list(identity, exp)) {
      ge <- group_effects(bg, summary = TRUE, transform = tr, prob = p)
      fe <- fixed_effects(bg, summary = TRUE, transform = tr, prob = p)
      expect_equal(unname(ge[1,]),
                   unname(expected_interval(tr(draws$theta_k[,1]), p)))
      expect_equal(unname(fe[1,]),
                   unname(expected_interval(tr(draws$beta[,1]), p)))
      te <- treatment_effect(bg, summary = TRUE, transform = tr, prob = p)
      expect_equal(unname(te$theta[1,]),
                   unname(expected_interval(tr(draws$theta), p)))
      expect_true(all(is.na(te$sigma[, -6])))
      expect_equal(unname(te$sigma[1,"width"]), p)
    }
    set.seed(130)
    de <- effect_draw(bg)
    set.seed(130)
    got <- effect_draw(bg, summary = TRUE, transform = exp, prob = p)
    expect_equal(unname(got[1,]),
                 unname(expected_interval(exp(de), p)))
    sel <- selection(bg, prob = p)
    expect_equal(unname(baggr:::interval_stats(sel, p)[1,]),
                 unname(expected_interval(draws$omega[,1], p)))
    cor <- mutau_cor(bg, summary = TRUE, prob = p)
    expect_equal(unname(cor[1,]),
                 unname(expected_interval(bg$fit[,2], p)))
    pool_draws <- heterogeneity(bg, summary = FALSE)[,1,1]
    expect_equal(unname(heterogeneity(bg, prob = p)[1,]),
                 unname(expected_interval(pool_draws, p)))
  }
})

test_that("printing selects named quantities even when summary order changes", {
  mock_interval_extract()
  original_group <- group_effects
  original_fixed <- fixed_effects
  original_mint <- mint
  original_pool <- pooling
  testthat::local_mocked_bindings(
    group_effects = function(...) {
      x <- testthat::with_mocked_bindings(original_group(...),
        mint = original_mint, fixed_effects = original_fixed, .package = "baggr")
      if(!"width" %in% colnames(x)) return(x)
      cols <- c("sd","upper","width","mean","lower","median")
      if(is.matrix(x)) x[,cols,drop=FALSE] else x[,cols,,drop=FALSE]
    },
    fixed_effects = function(...) {
      x <- testthat::with_mocked_bindings(original_fixed(...),
        mint = original_mint, .package = "baggr")
      if(!"width" %in% colnames(x)) return(x)
      cols <- c("median","upper","sd","width","lower","mean")
      if(is.matrix(x)) x[,cols,drop=FALSE] else x[,cols,,drop=FALSE]
    },
    mint = function(...) {
      x <- original_mint(...)
      if(is.null(dim(x))) c(unrelated = 123456, rev(x)) else
        cbind(unrelated = 123456, x[,rev(colnames(x)),drop = FALSE])
    },
    pooling = function(...) {
      x <- testthat::with_mocked_bindings(original_pool(...),
        mint = original_mint, .package = "baggr")
      if(length(dim(x)) == 3L) return(x)
      x[,rev(colnames(x)),drop=FALSE]
    }, .package = "baggr")
  for(k in c(1,2)) {
    bg <- interval_fixture(k)
    draws <- attr(bg$fit, "draws")
    for(p in c(.05,.5,.8,.95)) {
      for(exponent in c(FALSE, TRUE)) {
        tr <- if(exponent) exp else identity
        output <- capture.output(print(bg, prob = p, exponent = exponent,
                                       group = TRUE, digits = 7))
        tab <- t(vapply(seq_len(k), function(i)
          expected_interval(tr(draws$theta_k[,i]), p), numeric(6)))
        rownames(tab) <- bg$data$group
        tab <- tab[,colnames(tab) != "width",drop=FALSE]
        if(exponent) tab <- tab[,setdiff(colnames(tab), "sd"),drop = FALSE]
        colnames(tab)[match(c("lower","upper"), colnames(tab))] <-
          baggr:::interval_names(p)
        pd <- original_pool(bg, summary = FALSE)
        tab <- cbind(tab, pooling = apply(pd[,,1,drop = FALSE], 2, mean))
        expect_printed_table(output, tab)
        fixed_tab <- matrix(expected_interval(tr(draws$beta[,1]), p), 1,
                            dimnames = list("x", names(expected_interval(1:3,p))))
        fixed_tab <- fixed_tab[,colnames(fixed_tab) != "width",drop=FALSE]
        if(exponent) fixed_tab <- fixed_tab[,setdiff(colnames(fixed_tab), "sd"),drop = FALSE]
        colnames(fixed_tab)[match(c("lower","upper"), colnames(fixed_tab))] <-
          baggr:::interval_names(p)
        expect_printed_table(output, fixed_tab)
        pool <- baggr:::interval_stats(original_pool(bg, type = "total", prob = p)[1,], p)
        expected_line <- paste("Total pooling (1 - I^2) =", format(pool["mean"], digits=7,trim=TRUE),
                               "with", paste0(100*p,"%"), "interval",
                               format(pool["lower"],digits=7,trim=TRUE), "to",
                               format(pool["upper"],digits=7,trim=TRUE))
        expect_true(any(startsWith(output, expected_line)))
        expect_false(any(grepl("123456", output)))
        omega <- expected_interval(draws$omega[,1],p)
        expect_true(any(grepl(paste0("with ",100*p,"% interval ",
          format(omega["lower"],digits=7,trim=TRUE)," to ",
          format(omega["upper"],digits=7,trim=TRUE)), output, fixed=TRUE)))
      }
    }
  }
})

test_that("forest plots use means and requested endpoints in every row", {
  mock_interval_extract()
  bg <- interval_fixture()
  z <- attr(bg$fit, "draws")
  for(p in c(.05,.5,.8,.95)) {
    for(show in c("inputs", "posterior", "both", "covariates")) {
      fp <- forest_plot(bg, show = show, prob = p)
      expect_equal(unname(tail(fp$estimates[,1,1],1)), mean(z$theta))
      expect_equal(unname(tail(fp$estimates[,2,1],1)),
                   unname(quantile(z$theta,(1-p)/2)))
      expect_equal(unname(tail(fp$estimates[,3,1],1)),
                   unname(quantile(z$theta,(1+p)/2)))
      raw <- if(show %in% c("inputs","both")) bg$data$tau else
        if(show == "posterior") colMeans(z$theta_k) else mean(z$beta)
      lower <- if(show %in% c("inputs","both"))
        bg$data$tau - qnorm((1+p)/2)*bg$data$se else
        if(show == "posterior") apply(z$theta_k,2,quantile,probs=(1-p)/2) else
          quantile(z$beta,(1-p)/2)
      expect_equal(unname(fp$estimates[seq_along(raw)+1,1,1]), unname(raw))
      expect_equal(unname(fp$estimates[seq_along(raw)+1,2,1]), unname(lower))
    }
  }
})

test_that("comparison and baggr plots retain quantity and model identities", {
  mock_interval_extract()
  bg <- interval_fixture()
  draws <- attr(bg$fit,"draws")
  for(p in c(.05,.5,.8,.95)) {
    cmp <- baggr_compare(A = bg, B = bg, prob = p)
    expect_equal(cmp$covariates$lower, rep(unname(quantile(draws$beta,(1-p)/2)),2))
    expect_equal(cmp$covariates$median, rep(median(draws$beta),2))
    expect_equal(rownames(cmp$mean_trt), c("A","B"))
    out <- capture.output(print(cmp, digits = 7))
    expect_printed_table(out, baggr:::interval_table(cmp$mean_trt,p))
    gp <- plot(cmp, compare="groups", grid_parameters=FALSE, prob=p,
               order=FALSE, add_values=TRUE, values_digits=7)[[1]]
    expect_equal(sort(unique(gp$data$model)), c("A","B"))
    expected_group <- expected_interval(draws$theta_k[,1],p)
    row <- gp$data[gp$data$model=="A" & gp$data$group=="Study 1",]
    expect_equal(row$median, unname(expected_group["median"]))
    expect_equal(row$lower, unname(expected_group["lower"]))
    built <- ggplot2::ggplot_build(gp)
    expect_true(any(vapply(built$data, function(d)
      "label" %in% names(d) && any(grepl(as.character(signif(row$median,7)), d$label, fixed=TRUE)),
      logical(1))))
    for(values_outer in c(FALSE,TRUE)) {
      bp <- baggr_plot(bg, style="forest_plot", hyper=TRUE, order=FALSE,
                       prob=min(.25,p), prob_outer=p, values_outer=values_outer,
                       values_digits=7)
      row <- bp$data[bp$data$parameter=="Study 1",]
      expect_equal(row$m, median(draws$theta_k[,1]))
      expect_equal(row$ll, unname(quantile(draws$theta_k[,1],(1-p)/2)))
      expect_equal(row$hh, unname(quantile(draws$theta_k[,1],(1+p)/2)))
      built <- ggplot2::ggplot_build(bp)
      expect_equal(sort(built$data[[2]]$xmin), sort(bp$data$ll))
      expect_equal(sort(built$data[[3]]$xmin), sort(bp$data$l))
      text <- built$data[[4]]$label
      expected_text <- baggr:::paste3_formatter(row$m,
        if(values_outer) row$ll else row$l, if(values_outer) row$hh else row$h,7)
      expect_true(expected_text %in% text)
    }
    calls <- list()
    testthat::with_mocked_bindings({
      testthat::with_mocked_bindings(
        plot(cmp, grid_models=TRUE, prob=p),
        grid.arrange = function(grobs, ncol) invisible(grobs),
        .package = "gridExtra")
    }, baggr_plot = function(bg, prob, prob_outer, ...) {
      calls[[length(calls)+1]] <<- c(prob,prob_outer)
      ggplot2::ggplot()
    }, .package="baggr")
    expect_equal(calls, rep(list(c(min(.5,p),p)),2))
  }
})

test_that("pooling extremes and multiple effects keep names and dimensions", {
  mock_interval_extract()
  for(k in c(1,2)) {
    for(pool in c("none","full")) {
      bg <- interval_fixture(k, covariates = FALSE)
      bg$pooling <- pool
      bg$inputs$M <- 0
      for(p in c(0,.5,1)) {
        expect_output(print(bg, prob=p, digits=7), "Aggregate treatment effect")
        po <- pooling(bg, prob=p)
        expect_equal(as.numeric(po[,"mean"]), rep(if(pool=="full") 1 else 0,k))
      }
    }
  }
  bg <- interval_fixture(covariates=FALSE)
  bg$pooling <- "full"
  bg$inputs$M <- 0
  bg$effects <- c("first","second")
  bg$n_parameters <- 2
  draws <- attr(bg$fit,"draws")
  draws$theta <- cbind(first=draws$theta, second=draws$theta*2+1)
  attr(bg$fit,"draws") <- draws
  for(p in c(.5,.8)) {
    ge <- group_effects(bg, summary=TRUE, prob=p, transform=exp)
    expect_equal(dim(ge), c(2L,6L,2L))
    expect_equal(dimnames(ge)[[3]], bg$effects)
    expect_equal(ge[1,"mean","second"], mean(exp(draws$theta[,"second"])))
    output <- capture.output(print(bg, prob=p, digits=7))
    tab <- t(apply(draws$theta,2,expected_interval,prob=p))
    tab <- tab[,colnames(tab) != "width",drop=FALSE]
    colnames(tab)[match(c("lower","upper"),colnames(tab))] <- baggr:::interval_names(p)
    expect_printed_table(output,tab)
    plots <- baggr_plot(bg, style="forest_plot", prob=.25, prob_outer=p,
                        order=FALSE, hyper=TRUE)
    expect_length(plots,2)
    expect_equal(plots[[2]]$data$m, rep(median(draws$theta[,"second"]),3))
  }
})

test_that("aggregate and correlation printouts select the right quantities", {
  mock_interval_extract()
  bg <- interval_fixture(covariates=FALSE)
  bg$model <- "vivs"
  bg$inputs$M <- 0
  draws <- attr(bg$fit,"draws")
  draws$theta <- array(c(draws$theta*0,draws$theta), c(length(draws$theta),1,2))
  draws[["sigma[1,2]"]] <- draws$sigma
  attr(bg$fit,"draws") <- draws
  check_line <- function(output, label, vals, p) {
    line <- output[startsWith(output, label)]
    expect_length(line, 1)
    line <- substring(line, nchar(label) + 1)
    nums <- regmatches(line, gregexpr("[-+]?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?", line))[[1]]
    expect_equal(as.numeric(nums), unname(c(vals["mean"],100*p,vals["lower"],vals["upper"])),
                 tolerance = 1e-6)
  }
  for(p in c(.05,.5,.8,.95)) {
    set.seed(130)
    pred <- effect_draw(bg)
    set.seed(130)
    out <- capture.output(print(bg, prob=p, group=FALSE, fixed=FALSE,digits=7))
    check_line(out, "Hyper-SD (sigma) =", expected_interval(draws$sigma,p),p)
    check_line(out, "Posterior predictive effect =", expected_interval(pred,p),p)
    check_line(out, "Correlation of treatment effect and baseline =",
               expected_interval(bg$fit[,2],p),p)
    check_line(out, "Hypermean (theta) =", expected_interval(draws$theta[,1,2],p),p)
  }
})

test_that("comparison covariates retain the model name when other models omit them", {
  mock_interval_extract()
  bg <- interval_fixture()
  bg_none <- interval_fixture(covariates=FALSE)
  cmp <- suppressMessages(baggr_compare(no_covariates=bg_none, with_covariates=bg,
                                        prob=.5))
  expect_identical(cmp$covariates$model, "with_covariates")
  expect_equal(cmp$covariates$mean, mean(attr(bg$fit,"draws")$beta))
})

test_that("forest annotations follow the selected input or posterior quantity", {
  mock_interval_extract()
  bg <- interval_fixture()
  draws <- attr(bg$fit,"draws")
  for(show in c("inputs","posterior","both")) {
    for(which in c("inputs","posterior")) {
      fp <- forest_plot(bg,show=show,print=which,prob=.5,digits=7)
      means <- if(which=="inputs") bg$data$tau else colMeans(draws$theta_k)
      sds <- if(which=="inputs") bg$data$se else apply(draws$theta_k,2,sd)
      expect_equal(unlist(fp$labels[[2]][2:3]),
        paste0(format(means,digits=7)," (",format(sds,digits=7),")"))
      expect_equal(unlist(fp$labels[[1]][2:3]), bg$data$group)
    }
  }
  expect_error(forest_plot(bg,show="covariates",print="inputs"), "require")
})

test_that("forest plots draw means outside narrow central intervals", {
  mock_interval_extract()
  bg <- interval_fixture()
  theta <- attr(bg$fit, "draws")$theta
  drawn <- list()
  drawn_positions <- numeric()
  draw_ci <- function(lower_limit, estimate, upper_limit, ...) {
    drawn[[length(drawn) + 1]] <<- c(lower = lower_limit, mean = estimate,
                                    upper = upper_limit)
    drawn_positions <<- c(drawn_positions, grid::convertX(
      grid::unit(estimate, "native"), "npc", valueOnly = TRUE))
    forestplot::fpDrawNormalCI(lower_limit, estimate, upper_limit, ...)
  }
  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path)
  tryCatch({
    for(xlog in c(FALSE, TRUE)) {
      drawn <- list()
      drawn_positions <- numeric()
      fp <- forest_plot(bg, show = "posterior", prob = .05,
                        fn.ci_norm = draw_ci, boxsize = .15, xlog = xlog)
      expected <- c(lower = unname(quantile(theta, .475)), mean = mean(theta),
                    upper = unname(quantile(theta, .525)))
      if(xlog) expected <- log(expected)
      expect_gt(expected["mean"], expected["upper"])
      expect_equal(unname(fp$estimates[5,c("lower","mean","upper"),1]),
                   unname(expected))
      expect_gte(max(fp$xticks), expected["mean"])
      print(fp)
      expect_length(drawn, bg$n_groups + 1)
      expect_true(all(drawn_positions >= 0 & drawn_positions <= 1))
      expect_equal(unname(drawn[[length(drawn)]]), unname(expected))
    }
    expect_error(forest_plot(bg, prob = .05, xlog = TRUE, clip = c(-1, 2)),
                 "exponential form")
  }, finally = {
    grDevices::dev.off()
    unlink(plot_path)
  })
})

# Check the common contract at each public summary entry point.
expect_summary_contract <- function(x, rows, effects, prob) {
  stats <- c("mean", "median", "sd", "lower", "upper", "width")
  expect_type(x, "double")
  expect_identical(colnames(x), stats)
  if(length(effects) == 1) {
    expect_identical(dim(x), c(as.integer(rows), 6L))
    expect_identical(attr(x, "effect"), effects)
  } else {
    expect_identical(dim(x), c(as.integer(rows), 6L, length(effects)))
    expect_identical(dimnames(x)[[3]], effects)
  }
  for(i in seq_along(effects)) {
    mat <- baggr:::effect_slice(x, i)
    expect_equal(unname(mat[,"width"]), rep(prob, rows))
    expect_identical(dim(mat), c(as.integer(rows), 6L))
    expect_identical(attr(mat, "effect"), effects[i])
  }
}

test_that("single and multiple rows and effects follow the numeric contract", {
  mock_interval_extract()
  for(k in c(1, 2)) {
    bg <- interval_fixture(k)
    for(p in c(0, .05, .8, 1)) {
      for(f in list(group_effects, study_effects, random_effects, pooling, weights)) {
        x <- f(bg, summary = TRUE, prob = p)
        expect_summary_contract(x, k, bg$effects, p)
        expect_identical(rownames(x), bg$data$group)
      }
      x <- fixed_effects(bg, summary = TRUE, prob = p)
      expect_summary_contract(x, 1, bg$effects, p)
      expect_identical(rownames(x), "x")
      expect_summary_contract(selection(bg, prob = p), 1, "selection", p)
      expect_summary_contract(mutau_cor(bg, summary = TRUE, prob = p),
                              1, "correlation", p)
      expect_summary_contract(heterogeneity(bg, prob = p), 1, bg$effects, p)
      for(f in list(hypermean, hypersd, effect_draw, predict))
        expect_summary_contract(f(bg, summary = TRUE, prob = p), 1, bg$effects, p)
      newdata <- data.frame(x = c(-1, 2), row.names = c("new A", "new B"))
      for(nd in list(newdata[1,,drop = FALSE], newdata)) {
        set.seed(130)
        raw <- effect_draw(bg, newdata = nd)
        set.seed(130)
        x <- predict(bg, newdata = nd, summary = TRUE, prob = p)
        expect_summary_contract(x, nrow(nd), bg$effects, p)
        expect_identical(rownames(x), rownames(nd))
        raw <- array(raw, c(1, nrow(nd), 101))
        for(j in seq_len(nrow(nd)))
          expect_equal(unname(x[j,]), unname(expected_interval(raw[1,j,], p)))
      }
    }
  }
  bg <- interval_fixture(covariates = FALSE)
  bg$effects <- c("first", "second")
  bg$n_parameters <- 2
  draws <- attr(bg$fit, "draws")
  draws$theta <- cbind(draws$theta, draws$theta*2)
  draws$sigma <- cbind(draws$sigma, draws$sigma*2)
  draws$theta_k <- array(c(draws$theta_k, draws$theta_k*2), c(101, 2, 2))
  attr(bg$fit, "draws") <- draws
  for(pool in c("partial", "full", "none")) {
    bg$pooling <- pool
    for(p in c(0, .8, 1)) {
      te <- treatment_effect(bg, summary = TRUE, prob = p, message = FALSE)
      for(x in te) expect_summary_contract(x, 1, bg$effects, p)
      expect_summary_contract(group_effects(bg, summary = TRUE, prob = p),
                              2, bg$effects, p)
      expect_summary_contract(pooling(bg, prob = p), 2, bg$effects, p)
      expect_summary_contract(heterogeneity(bg, prob = p), 1, bg$effects, p)
      expect_summary_contract(effect_draw(bg, summary = TRUE, prob = p,
                                         message = FALSE), 1, bg$effects, p)
      expect_summary_contract(fixed_effects(bg, summary = TRUE, prob = p),
                              0, bg$effects, p)
      if(pool %in% c("none", "full")) {
        cor <- mutau_cor(bg, summary = TRUE, prob = p)
        expect_summary_contract(cor, 1, "correlation", p)
        expect_true(all(is.na(cor[,-6])))
      }
      if(pool == "full") expect_equal(unname(te$sigma[,"sd",]), c(0, 0))
      if(pool == "none") expect_true(all(is.na(te$theta[,-6,,drop = FALSE])))
      transformed <- treatment_effect(bg, summary = TRUE, transform = exp,
                                       prob = p, message = FALSE)$sigma
      expect_summary_contract(transformed, 1, bg$effects, p)
      expect_true(all(is.na(transformed[,-6,,drop = FALSE])))
    }
    cmp <- suppressMessages(baggr_compare(A = bg, B = bg, prob = .8))
    for(x in cmp[c("mean_trt", "sd_trt", "posteriorpd_trt")]) {
      expect_summary_contract(x, 2, bg$effects, .8)
      expect_identical(rownames(x), c("A", "B"))
    }
    expect_no_error(capture.output(print(cmp)))
  }
  bg$model <- "quantiles"
  bg$pooling <- "partial"
  expect_warning(x <- heterogeneity(bg, prob = .8), "not implemented")
  expect_summary_contract(x, 1, bg$effects, .8)
  expect_true(all(is.na(x[,-6,,drop = FALSE])))
})

test_that("percentage mode changes only endpoint names and drops width", {
  mock_interval_extract()
  for(n_effects in c(1, 2)) {
    bg <- interval_fixture()
    if(n_effects == 2) {
      bg$effects <- c("first", "second")
      bg$n_parameters <- 2
      bg$pooling <- "full"
      draws <- attr(bg$fit, "draws")
      draws$theta <- cbind(draws$theta, draws$theta*2)
      attr(bg$fit, "draws") <- draws
    }
    for(p in c(0, .05, .8, 1)) {
      for(f in list(group_effects, study_effects, random_effects,
                    treatment_effect, hypermean, hypersd)) {
        a <- f(bg, summary = TRUE, prob = p)
        b <- f(bg, summary = TRUE, prob = p, rename_int = TRUE)
        if(is.list(a)) {
          a <- a[c("theta", "sigma")]
          b <- b[c("theta", "sigma")]
        } else { a <- list(a); b <- list(b) }
        for(j in seq_along(a)) {
          expect_identical(colnames(b[[j]]),
            c("mean", "median", "sd", baggr:::interval_names(p)))
          for(i in seq_len(n_effects)) {
            plain <- baggr:::effect_slice(a[[j]], i)
            percent <- baggr:::effect_slice(b[[j]], i)
            expect_identical(as.numeric(percent),
                             as.numeric(plain[, -6, drop = FALSE]))
            expect_identical(attr(percent, "effect"), attr(plain, "effect"))
          }
        }
      }
    }
  }
  expect_equal(mint(rep(2, 10), sd = TRUE)["sd"], c(sd = 0))
  expect_equal(mint(NA_real_, median = TRUE, sd = TRUE),
    c(mean = NA_real_, median = NA_real_, sd = NA_real_, lower = NA_real_,
      upper = NA_real_, width = .95))
  expect_identical(mint(1:9, digits = 2)["width"], c(width = .95))
  expect_identical(dim(mint(matrix(1:9, ncol = 1))), c(1L, 4L))
})

test_that("LOO plots select named medians and bounds after reordering", {
  mock_interval_extract()
  original <- treatment_effect
  testthat::local_mocked_bindings(treatment_effect = function(...) {
    x <- original(...)
    x$theta <- x$theta[,c("upper", "sd", "width", "lower", "mean", "median"),
                   drop = FALSE]
    x
  }, .package = "baggr")
  bg <- interval_fixture()
  loo <- structure(list(full_model = bg, models = list(bg, bg),
                        pointwise = c(-1, -2)), class = "baggr_cv")
  gp <- plot(loo, add_values = FALSE)
  vals <- expected_interval(attr(bg$fit, "draws")$theta, .95)
  got <- gp$data[gp$data$model == "LOOCV estimate",]
  for(stat in c("median", "lower", "upper"))
    expect_equal(got[[stat]], rep(unname(vals[stat]), 2))
  expect_no_error(ggplot2::ggplot_build(gp))
})

test_that("single-study VIVS extraction preserves its summary row", {
  mock_interval_extract()
  bg <- interval_fixture(1, covariates = FALSE)
  bg$model <- "vivs"
  draws <- attr(bg$fit, "draws")
  draws$theta_k <- array(c(draws$theta_k*0, draws$theta_k), c(101, 1, 2, 1))
  draws$theta <- array(c(draws$theta*0, draws$theta), c(101, 1, 2))
  draws[["sigma[1,2]"]] <- array(draws$sigma, 101)
  attr(bg$fit, "draws") <- draws
  got <- group_effects(bg, summary = TRUE, prob = .8)
  expect_summary_contract(got, 1, bg$effects, .8)
  expect_identical(rownames(got), "Study 1")
  expect_equal(unname(got[1,]), unname(expected_interval(draws$theta_k[,1,2,1], .8)))
  te <- treatment_effect(bg, summary = TRUE, prob = .8)
  for(x in te) expect_summary_contract(x, 1, bg$effects, .8)
  expect_equal(unname(te$sigma[1,]),
               unname(expected_interval(draws$sigma, .8)))
  expect_identical(dim(treatment_effect(bg)$sigma), 101L)

})

test_that("fixed, quantile and LOO consumers slice multiple effects by name", {
  mock_interval_extract()
  bg <- interval_fixture()
  bg$effects <- c("first", "second")
  bg$n_parameters <- 2
  draws <- attr(bg$fit, "draws")
  draws$theta <- cbind(draws$theta, draws$theta*2)
  draws$sigma <- cbind(draws$sigma, draws$sigma*2)
  draws$theta_k <- array(c(draws$theta_k, draws$theta_k*2), c(101, 2, 2))
  attr(bg$fit, "draws") <- draws
  fe <- fixed_effects(bg, summary = TRUE, prob = .8)
  expect_summary_contract(fe, 1, bg$effects, .8)
  expect_identical(rownames(fe), "x")
  expect_identical(baggr:::effect_slice(fe, "second"), baggr:::effect_slice(fe, 2))
  expect_equal(unname(fe[1,,2]), unname(expected_interval(draws$beta[,1], .8)))
  bg$summary_data <- rbind(data.frame(bg$data, treatment = "B"),
                           data.frame(bg$data, treatment = "C"))
  attr(bg$inputs, "treatment_levels") <- c("B", "C")
  original_group <- group_effects
  original_treatment <- treatment_effect
  testthat::local_mocked_bindings(
    group_effects = function(...) {
      x <- original_group(...)
      if("width" %in% colnames(x)) x[,rev(colnames(x)),,drop = FALSE] else x
    },
    treatment_effect = function(...) {
      x <- original_treatment(...)
      if("width" %in% colnames(x$theta))
        x$theta <- x$theta[,rev(colnames(x$theta)),,drop = FALSE]
      x
    }, .package = "baggr")
  loo <- structure(list(full_model = bg, models = list(bg, bg),
                        pointwise = c(-1, -2)), class = "baggr_cv")
  gp <- plot(loo, add_values = FALSE)
  for(i in seq_along(bg$effects)) {
    got <- gp$data[gp$data$model == "LOOCV estimate" &
                   gp$data$parameter == bg$effects[i],]
    expect_equal(got$median, rep(median(draws$theta[,i]), 2))
  }
  expect_no_error(ggplot2::ggplot_build(gp))
  bg$covariates <- character()
  output <- capture.output(print(bg, digits = 7))
  tab <- t(vapply(seq_len(2), function(j)
    expected_interval(draws$theta_k[,j,2], .95), numeric(6)))
  rownames(tab) <- bg$data$group
  tab <- baggr:::interval_table(tab, .95)
  tab <- cbind(tab, pooling = pooling(bg)[,"mean",2])
  expect_printed_table(output, tab)
  cmp <- baggr_compare(A = bg, B = bg)
  expect_equal(unname(cmp$mean_trt[,"mean",2]), rep(mean(draws$theta[,2]), 2))
  gp_cmp <- plot(cmp, grid_parameters = FALSE, order = FALSE)[[2]]
  expect_equal(gp_cmp$data$median[gp_cmp$data$group == "Study 1"],
               rep(median(draws$theta_k[,1,2]), 2))
  expect_no_error(ggplot2::ggplot_build(gp_cmp))
  ge <- group_effects(bg, summary = TRUE)
  te <- treatment_effect(bg, summary = TRUE)
  bg$model <- "quantiles"
  bg$quantiles <- c(.25, .75)
  plots <- testthat::with_mocked_bindings(plot_quantiles(bg),
    group_effects = function(...) ge, treatment_effect = function(...) te,
    .package = "baggr")
  expect_equal(plots[[1]]$data$mean, unname(colMeans(draws$theta)))
  expect_equal(plots[[2]]$data$mean,
               c(colMeans(draws$theta_k[,,1]), colMeans(draws$theta_k[,,2])))
  expect_no_error(lapply(plots, ggplot2::ggplot_build))
})
