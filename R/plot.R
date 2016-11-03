plotMotifCount = function(object, percentage = FALSE, ...) {
  x = getMotifCount(object, percentage)
  mp = barplot(x, las = 1, axes = FALSE, axisnames = FALSE, ...)
  axis(2, las = 1)
  mtext(1:20, side = 1, at = mp, cex = 0.8)
  box()
}

plotCounts = function(x, cut) {
  plot(x$Counts, ylim = c(0, 100), xlab = "Architectures", ylab = "Counts/Percentage", axes = FALSE, type = "l")
  points(x$Counts, col = "black", pch = 21, bg = "gray")
  lines(x$Percentage, col = "darkblue")
  lines(x$Cumulative, col = "darkred")
  points(x$Percentage, col = "darkblue", pch = 21, bg = "steelblue")
  points(x$Cumulative, col = "darkred", pch = 21, bg = "orange")
  if(!missing(cut))
    abline(v = cut, lty = "dotted")
  axis(2, las = 1)
  legend("right", c("Counts", "Percentage", "Cumulative"), pch = 21, col = c("black", "darkblue", "darkred"), pt.bg = c("gray", "steelblue", "orange"), bty = "n")
  box()
}