obs_pred <- function(bg, interval = .8) {
  # rbind(
  #   as.data.frame(group_effects(bg, interval = interval), type = "model"),
  #
  #
  #
  # ggplot(plotdata, aes(x = observed, y = median, ymin = lb, ymax = ub)) +
  #   geom_hline(yintercept = tot_avg, color = "lightpink", size = 0.75) +
  #   geom_abline(intercept = 0, slope = 1, color = "skyblue") +
  #   geom_linerange(color = "gray60", size = 0.75) +
  #   geom_point(size = 2.5, shape = 21, fill = "gray30", color = "white", stroke = 0.2) +
  #   facet_grid(. ~ model) +
  #   coord_fixed() +
  #   scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) +
  #   labs(x = "Observed Hits / AB", y = "Predicted chance of hit") +
  #   ggtitle("Posterior Medians and 80% Intervals")
}
