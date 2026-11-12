# Two-sample Mendelian randomization: BMI -> coronary heart disease
# Exposure: ieu-a-2 (BMI, GIANT); Outcome: ieu-a-7 (CHD, CARDIoGRAMplusC4D)
# Requires an OpenGWAS token in ~/.Renviron as OPENGWAS_JWT (never commit it)

library(TwoSampleMR)

exp <- extract_instruments(outcomes = "ieu-a-2")
out <- extract_outcome_data(snps = exp$SNP, outcomes = "ieu-a-7")
dat <- harmonise_data(exp, out)

res <- mr(dat)
write.csv(res, "results/mr_results.csv", row.names = FALSE)
write.csv(mr_heterogeneity(dat), "results/heterogeneity.csv", row.names = FALSE)
write.csv(mr_pleiotropy_test(dat), "results/pleiotropy.csv", row.names = FALSE)

ggplot2::ggsave("results/scatter.png", mr_scatter_plot(res, dat)[[1]], width = 7, height = 7)
ggplot2::ggsave("results/forest.png", mr_forest_plot(mr_singlesnp(dat))[[1]], width = 7, height = 8)
ggplot2::ggsave("results/leaveoneout.png", mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]], width = 7, height = 8)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
