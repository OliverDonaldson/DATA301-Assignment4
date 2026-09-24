suppressMessages({library(fpp2); library(ggplot2)})
p <- function(g, f, w=9, h=4.5) ggsave(f, g, width=w, height=h, dpi=110)
p(autoplot(ausbeer) + labs(title="Australian quarterly beer production, 1956 Q1 – 2010 Q2",
    y="Megalitres", x="Year"), "e2_time.png")
p(ggseasonplot(ausbeer, year.labels=TRUE, year.labels.left=TRUE) +
    labs(title="Seasonal plot", y="Megalitres"), "e2_season.png", 9, 5.5)
p(ggseasonplot(ausbeer, polar=TRUE) + labs(title="Polar seasonal plot", y="Megalitres"), "e2_polar.png", 7, 5.5)
p(ggsubseriesplot(ausbeer) + labs(title="Seasonal subseries plot", y="Megalitres"), "e2_subseries.png")
p(ggAcf(ausbeer, lag.max=24) + labs(title="ACF of ausbeer"), "e2_acf.png", 9, 3.6)
# complete years only, one row per year
full <- window(ausbeer, end=c(2009,4))
m <- matrix(full, ncol=4, byrow=TRUE)
d <- data.frame(year=1956:2009, level=rowMeans(m), amp=m[,4]-m[,2])
d$decade <- paste0(d$year %/% 10 * 10, "s")
p(autoplot(aggregate(full, FUN=sum)) +
    labs(title="Annual totals (complete years 1956-2009)", subtitle="Peak 1981; the seasonal swing is removed by aggregation",
         x="Year", y="Megalitres"), "e2_annual.png", 8, 3.8)
# does the swing grow with the level? (a plot only: Exercise 2 fits no model)
p(ggplot(d, aes(level, amp, colour=decade)) + geom_point(size=1.8) +
    labs(title="Seasonal amplitude (Q4 - Q2) against annual level",
         subtitle="Rises with the level to the 1980s; the 2000s drop at almost the 1990s level",
         x="Annual mean (ML)", y="Q4 - Q2 (ML)", colour=NULL), "e2_amp.png", 7, 4.2)
# the Q8 visual argument: raw vs log
b <- rbind(data.frame(t=as.numeric(time(ausbeer)), y=as.numeric(ausbeer), s="Original (ML)"),
           data.frame(t=as.numeric(time(ausbeer)), y=log(as.numeric(ausbeer)), s="Log scale"))
b$s <- factor(b$s, levels=c("Original (ML)","Log scale"))
p(ggplot(b, aes(t, y)) + geom_line(linewidth=.35) + facet_wrap(~s, ncol=1, scales="free_y") +
    labs(title="Would a log transformation help?",
         subtitle="Log: the swing shrinks over time - early growth over-corrected, the 2000s contraction remains",
         x="Year", y=NULL), "e2_logcheck.png", 9, 5)
# decade swing on several Box-Cox scales: no power makes it constant
for (l in c(1, .5, BoxCox.lambda(ausbeer), 0)) {
  w <- matrix(BoxCox(full, l), ncol=4, byrow=TRUE); s <- tapply(w[,4]-w[,2], d$year %/% 10, mean)
  cat(sprintf("lambda=%.2f  largest/smallest decade swing=%.2f\n", l, max(s)/min(s)))
}
