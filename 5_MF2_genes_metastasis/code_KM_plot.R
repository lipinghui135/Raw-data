library(survival)
library(survminer)

df = read.table("./CAMK2N1_COLEC11_survival_data.txt", row.names=1, header=T, sep="\t")
fit = survfit(Surv(time, status)~CAMK2N1_type, data=df)
p  =ggsurvplot(fit,
	data = df,
	conf.int = TRUE,
	pval = TRUE,
	surv.median.line = "hv",
	risk.table = TRUE,
	xlab = "Time(years)",
	palette = c("hotpink", "chartuse"),
	legend.labs = c("Higher Expression", "Lower Expression"),
	legend.title = "Risk",
	break.x.by = 3,
	pval.size = 7,
	font.legend = 12
)
ggsave("CAMK2N1_survial_plot.pdf", p, width=6, height=6)

### 
