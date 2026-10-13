# Exact output checks use asymmetric draws, without fitting a model.
interval_fixture <- function(n_groups = 2, covariates = TRUE) {
  z <- seq(-1, 2, length.out = 101)^3 / 4
  fit <- cbind(1, seq(-0.8, 0.6, length.out = 101))
  attr(fit, "draws") <- list(mu = z + 2, tau = exp(z)/5,
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
  c(mean = mean(x), sd = sd(x), median = median(x),
    lci = unname(quantile(x, (1 - prob)/2)),
    uci = unname(quantile(x, (1 + prob)/2)))
}

expect_printed_table <- function(output, tab) {
  expect_true(all(capture.output(print(tab, digits = 7)) %in% output))
}

test_that("prob means contained mass and legacy arguments remain compatible", {
  mock_interval_extract()
  bg <- interval_fixture()
  z <- attr(bg$fit, "draws")$mu
  for(p in c(0, .05, .5, .8, .95, 1)) {
    got <- mint(z, prob = p, median = TRUE, sd = TRUE)
    expect_equal(unname(got), unname(c(quantile(z, (1-p)/2), mean(z),
                                      quantile(z, (1+p)/2), median(z), sd(z))))
    expect_identical(mint(z, p), mint(z, int = p))
    for(f in list(group_effects, fixed_effects, treatment_effect,
                  hypermean, hypersd, selection, mutau_cor)) {
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
      expect_equal(unname(ge[1,c("mean","sd","median","lci","uci"),1]),
                   unname(expected_interval(tr(draws$theta_k[,1]), p)))
      expect_equal(unname(fe[1,c("mean","sd","median","lci","uci"),1]),
                   unname(expected_interval(tr(draws$beta[,1]), p)))
      te <- treatment_effect(bg, summary = TRUE, transform = tr, prob = p)
      expect_equal(unname(baggr:::interval_stats(te$tau, p)),
                   unname(expected_interval(tr(draws$mu), p)))
      expect_true(all(is.na(te$sigma_tau)))
    }
    set.seed(130)
    de <- effect_draw(bg)
    set.seed(130)
    got <- effect_draw(bg, summary = TRUE, transform = exp, prob = p)
    expect_equal(unname(baggr:::interval_stats(got, p)),
                 unname(expected_interval(exp(de), p)))
    sel <- selection(bg, prob = p)
    expect_equal(unname(baggr:::interval_stats(sel, p)[1,]),
                 unname(expected_interval(draws$omega[,1], p)))
    cor <- mutau_cor(bg, summary = TRUE, prob = p)
    expect_equal(unname(baggr:::interval_stats(cor, p)),
                 unname(expected_interval(bg$fit[,2], p)))
    pool_draws <- heterogeneity(bg, summary = FALSE)[,1,1]
    expect_equal(unname(heterogeneity(bg, prob = p)[,1,1]),
                 unname(c(quantile(pool_draws, (1-p)/2), mean(pool_draws),
                          quantile(pool_draws, (1+p)/2))))
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
      if(length(dim(x)) != 3L) return(x)
      x[,c("sd","uci","mean","lci","median"),,drop = FALSE]
    },
    fixed_effects = function(...) {
      x <- testthat::with_mocked_bindings(original_fixed(...),
        mint = original_mint, .package = "baggr")
      if(length(dim(x)) != 3L) return(x)
      x[,c("median","uci","sd","lci","mean"),,drop = FALSE]
    },
    mint = function(...) {
      x <- original_mint(...)
      if(is.null(dim(x))) c(unrelated = 123456, rev(x)) else
        cbind(unrelated = 123456, x[,rev(colnames(x)),drop = FALSE])
    },
    pooling = function(...) {
      x <- testthat::with_mocked_bindings(original_pool(...),
        mint = original_mint, .package = "baggr")
      x[rev(seq_len(dim(x)[1])),,,drop = FALSE]
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
          expected_interval(tr(draws$theta_k[,i]), p), numeric(5)))
        rownames(tab) <- bg$data$group
        if(exponent) tab <- tab[,setdiff(colnames(tab), "sd"),drop = FALSE]
        colnames(tab)[match(c("lci","uci"), colnames(tab))] <-
          baggr:::interval_names(p)
        pd <- original_pool(bg, summary = FALSE)
        tab <- cbind(tab, pooling = apply(pd[,,1,drop = FALSE], 2, mean))
        expect_printed_table(output, tab)
        fixed_tab <- matrix(expected_interval(tr(draws$beta[,1]), p), 1,
                            dimnames = list("x", names(expected_interval(1:3,p))))
        if(exponent) fixed_tab <- fixed_tab[,setdiff(colnames(fixed_tab), "sd"),drop = FALSE]
        colnames(fixed_tab)[match(c("lci","uci"), colnames(fixed_tab))] <-
          baggr:::interval_names(p)
        expect_printed_table(output, fixed_tab)
        pool <- baggr:::interval_stats(original_pool(bg, type = "total", prob = p)[,,1], p)
        expected_line <- paste("Total pooling (1 - I^2) =", format(pool["mean"], digits=7,trim=TRUE),
                               "with", paste0(100*p,"%"), "interval",
                               format(pool["lci"],digits=7,trim=TRUE), "to",
                               format(pool["uci"],digits=7,trim=TRUE))
        expect_true(any(startsWith(output, expected_line)))
        expect_false(any(grepl("123456", output)))
        omega <- expected_interval(draws$omega[,1],p)
        expect_true(any(grepl(paste0("with ",100*p,"% interval ",
          format(omega["lci"],digits=7,trim=TRUE)," to ",
          format(omega["uci"],digits=7,trim=TRUE)), output, fixed=TRUE)))
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
      expect_equal(unname(tail(fp$estimates[,1,1],1)), mean(z$mu))
      expect_equal(unname(tail(fp$estimates[,2,1],1)),
                   unname(quantile(z$mu,(1-p)/2)))
      expect_equal(unname(tail(fp$estimates[,3,1],1)),
                   unname(quantile(z$mu,(1+p)/2)))
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
    expect_equal(cmp$covariates$lci, rep(unname(quantile(draws$beta,(1-p)/2)),2))
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
    expect_equal(row$lci, unname(expected_group["lci"]))
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
        expect_equal(as.numeric(po["mean",,1]), rep(if(pool=="full") 1 else 0,k))
      }
    }
  }
  bg <- interval_fixture(covariates=FALSE)
  bg$pooling <- "full"
  bg$inputs$M <- 0
  bg$effects <- c("first","second")
  bg$n_parameters <- 2
  draws <- attr(bg$fit,"draws")
  draws$mu <- cbind(first=draws$mu, second=draws$mu*2+1)
  attr(bg$fit,"draws") <- draws
  for(p in c(.5,.8)) {
    ge <- group_effects(bg, summary=TRUE, prob=p, transform=exp)
    expect_equal(dim(ge), c(2L,5L,2L))
    expect_equal(dimnames(ge)[[3]], bg$effects)
    expect_equal(ge[1,"mean","second"], mean(exp(draws$mu[,"second"])))
    output <- capture.output(print(bg, prob=p, digits=7))
    tab <- t(apply(draws$mu,2,expected_interval,prob=p))
    colnames(tab)[match(c("lci","uci"),colnames(tab))] <- baggr:::interval_names(p)
    expect_printed_table(output,tab)
    plots <- baggr_plot(bg, style="forest_plot", prob=.25, prob_outer=p,
                        order=FALSE, hyper=TRUE)
    expect_length(plots,2)
    expect_equal(plots[[2]]$data$m, rep(median(draws$mu[,"second"]),3))
  }
})

test_that("aggregate and correlation printouts select the right quantities", {
  mock_interval_extract()
  bg <- interval_fixture(covariates=FALSE)
  bg$model <- "mutau"
  bg$inputs$M <- 0
  draws <- attr(bg$fit,"draws")
  draws$mu <- array(c(draws$mu*0,draws$mu), c(length(draws$mu),1,2))
  draws[["hypersd[1,2]"]] <- draws$tau
  attr(bg$fit,"draws") <- draws
  check_line <- function(output, label, vals, p) {
    line <- output[startsWith(output, label)]
    expect_length(line, 1)
    line <- substring(line, nchar(label) + 1)
    nums <- regmatches(line, gregexpr("[-+]?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?", line))[[1]]
    expect_equal(as.numeric(nums), unname(c(vals["mean"],100*p,vals["lci"],vals["uci"])),
                 tolerance = 1e-6)
  }
  for(p in c(.05,.5,.8,.95)) {
    set.seed(130)
    pred <- effect_draw(bg)
    set.seed(130)
    out <- capture.output(print(bg, prob=p, group=FALSE, fixed=FALSE,digits=7))
    check_line(out, "Hyper-SD (sigma_tau) =", expected_interval(draws$tau,p),p)
    check_line(out, "Posterior predictive effect =", expected_interval(pred,p),p)
    check_line(out, "Correlation of treatment effect and baseline =",
               expected_interval(bg$fit[,2],p),p)
    check_line(out, "Hypermean (tau) =", expected_interval(draws$mu[,1,2],p),p)
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
