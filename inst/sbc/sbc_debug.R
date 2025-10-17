job_bad1 <- 8130
job_bad2 <- 1531

## used for debugging
job1 <- testJob(5)
job1

names(job1)

job1$rank


job1_def <- makeJob(job_bad1)
job1_def$pars

res1 <- fit_exnex(base_data_final, job1_def, job1_def$instance, save_fit = TRUE)
res2 <- fit_exnex(base_data, job1_def, job1_def$instance, save_fit = TRUE)
res3 <- fit_exnex(base_data_final, job1_def, job1_def$instance, save_fit = TRUE)

res1

res2

res3

w1 <- learn_warmup_info(res1$fit$standata, res1$fit$stanfit)
w2 <- learn_warmup_info(res2$fit$standata, res2$fit$stanfit)
w3 <- learn_warmup_info(res3$fit$standata, res3$fit$stanfit)

w1$stepsize
w2$stepsize
w3$stepsize
w1$inv_metric
w2$inv_metric
w3$inv_metric


base_data_final2 <- base_data
base_data_final2$models <- modifyList(base_data_final2$models, warmup_info_by_model)


fit <- res1$fit
learn_warmup_info(fit$standata, fit$stanfit)
draw <- extract_draw(rstan::extract(fit$stanfit)[1:8], 1)
names(draw)
restore_draw_dims(fit$standata, draw)
draw$tau_log_beta_raw
draw

job1_def <- makeJob(5)
job1_def$pars

res1 <- fit_exnex(base_data, job1_def, job1_def$instance, save_fit = FALSE)

length(res1$draws)

ini <- base_data_final[[1]]$log2bayes_EXNEX$warmup_info[[50]]$draw
stan_rdump(names(ini), "init.data.R", envir = list2env(ini))

res1

res1$fit
fit_sum <- rstan::summary(res1$fit$stanfit)$summary

fit_sum[apply(sapply(params, grepl, x = rownames(fit_sum)), 1, any), c("n_eff", "Rhat")]

post_lp <- as.array(res1$fit$stanfit, pars = "lp__")
class(post_lp)

dim(post_lp)
m <- monitor(post_lp)

names(m)
class(m)
rownames(m)
m["lp__", "Bulk_ESS"]
as.numeric(m[1, c("Bulk_ESS", "Tail_ESS")])

job2_def <- makeJob(5)
job2_def$pars

res2 <- fit_exnex(base_data_final, job2_def, job2_def$instance, save_fit = TRUE)

print(res1$fit$stanfit, pars = "lp__")
print(res2$fit$stanfit, pars = "lp__")


dim(job1_def$instance$draw$draw_beta)
job1_def$instance$draw$draw_beta[1, 1, , ]

summary(job1$fit$stanfit)

names(job1)

fit <- job1$fit

job2 <- testJob(6)
job3 <- testJob(11)

job <- makeJob(1)

attributes(job)
names(job)

options(mc.cores = 2)

out <- fit_exnex(data = base_data, job = job, instance = job$instance)
