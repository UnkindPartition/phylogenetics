library(tidyverse)

d <- read_csv("trace.csv") %>%
  group_by(method) %>%
  arrange(step) %>%
  mutate(cumtime = cumsum(time)) %>%
  ungroup

info <- read_csv("info.csv") %>%
  pivot_wider(names_from="var", values_from="value")

lik <- read_csv("likprofile.csv")

ggplot(d, aes(step, ll, color=method)) +
  labs(x="Step", y="Log likelihood") +
  geom_hline(aes(yintercept = info$true_ll)) +
  ylim(max(info$true_ll-10,min(d$ll)),NA) +
  geom_line()
ggsave("loglik.png", width=12)

ggplot(d, aes(step, rate, color=method)) +
  labs(x="Step", y="Learning rate") +
  expand_limits(y=0) +
  geom_line()
ggsave("rate.png", width=12)

ggplot(d, aes(step, mse, color=method)) +
  labs(x="Step", y="Mean squared error") +
  expand_limits(y=0) +
  geom_line()
ggsave("mse.png", width=12)

ggplot(d, aes(step, cumtime, color=method)) +
  labs(x="Step", y="Accumulated time") +
  geom_line()
ggsave("time.png", width=12)

ggplot(d, aes(cumtime, ll, color=method)) +
  labs(x="Time", y="Log likelihood") +
  geom_hline(aes(yintercept = info$true_ll)) +
  ylim(max(info$true_ll-10,min(d$ll)),NA) +
  geom_line()
ggsave("loglik_by_time.png", width=12)

ggplot(lik, aes(alpha, ll, color=method)) +
  labs(x="opt -> true", y="Log likelihood") +
  geom_line()
ggsave("likprofile.png", width=12)
