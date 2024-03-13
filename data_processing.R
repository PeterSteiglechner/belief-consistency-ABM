library(tidyverse)
library(corclass)
library(igraph)

raw <- read_csv('/home/peter/PhD_ZMT/Project_coherentagents/coherentagents/ESS9e03_1.csv') # Use this to select other context variables later
dffull <-  raw %>% 
  # Filtering the cases -- cases with missing values on believes variables deleted:
  filter(freehms <= 5, gincdif <= 5, impcntr <= 4, 
         lrscale <= 10, euftf <= 10)%>%
  mutate(
    # Scaling data to values within [-1, 1]. 
    # v' = -1 + 2 * (v-m)/(M-m), where m is the lowest, M the highest possible answer
    across(c(freehms, gincdif), ~ -1 + 2 * (.x - 1)/(5-1)), 
    across(c(lrscale, euftf), ~ -1 + 2 * (.x - 0)/(10-0)), 
    impcntr = -1 + 2* (impcntr - 1)/(4-1),
    # Flipping some scales: some questions are asked in a "negative" sense: 
    # e.g. -1 --> "more immigrants" and 1 --> "less immigrants"
    across(c(freehms, gincdif, impcntr), ~ .x * -1)
  )
attitudenames <- c("freehms", "gincdif", "lrscale", "impcntr", "euftf")
cntrynames <- dffull$cntry |> unique()

write_ABM_data <- function(country, attitudenames) {
  if (!dir.exists(paste0("ns_",country))) {
    dir.create(paste0("ns_",country))
  } 
  countrydf <- dffull |> select(idno, cntry, all_of(attitudenames)) |> filter(cntry == country)
  cm <- countrydf |> select(all_of(attitudenames)) |> cca(filter.significance = TRUE, filter.value = 1)
  
  1:length(cm$modules) |> 
    map(function(i) cm$modules[[i]]$dtf |> as_tibble() |>
          mutate(group = i)) |>
    reduce(bind_rows) |>
    mutate(idno = sample(1:nrow(countrydf), nrow(countrydf), replace = F)) |>  
    arrange(idno) |> 
             # Note: Our the most advanced ABM supposes that the first column of data is 'idno' which is not read in agents, so we "create" here 'idno' containing just unuseful number -1. 
             # Note JL: Added random resample of idnos to avoid ordered network
    relocate(idno, .before = freehms) |>
    write_csv(paste0("ns_",country,"/itemsCCA1.csv"))
  1:length(cm$modules) |> 
    map(function(i) cm$modules[[i]]$cormat |> as_tibble() |> 
          mutate(item = attitudenames, group = i)) |> 
    reduce(bind_rows) |> 
    write_csv(paste0("ns_",country,"/correlationsCCA1.csv"))
}  


write_ABM_data("DE",attitudenames)
## The following writes out groups for all countries
# cntrynames |> map(\(x) write_ABM_data(x,attitudenames))


## Eigenvector Analysis

R <- read_csv("ns_DE/correlationsCCA3.csv")
X <- read_csv("ns_DE/itemsCCA3.csv")
X7 <- read_csv("ns_DE/itemsCCA7.csv")

countrydf <- dffull |> select(idno, cntry, all_of(attitudenames), prtvede2) |> filter(cntry == "DE")
cm <- countrydf |> select(all_of(attitudenames)) |> cca(filter.significance = FALSE)
cm7 <- countrydf |> select(all_of(attitudenames)) |> cca(filter.significance = TRUE)


countrydf$group3 = cm$membership

countrydf$group7 = cm7$membership

a = as.data.frame(table(countrydf$group7, countrydf$group3))
names(a) <- c("group7", "group3", "Count")
a |> ggplot(aes(x=group7, y=group3, size=Count, color=group3) ) + geom_point(alpha=0.4) +  
  scale_size_continuous(range = c(1, 30)) +  # Adjust the size range as needed
  geom_text(size = 3, vjust = -0.5) +  # Adjust text size and position as needed
  labs(title = "Grouping",
       x = "7-CCA grouping",
       y = "3-CCA grouping",
       size = "Count")

countrydf |> filter(group==1) |> select(prtvede2) |> table() / 1040
countrydf |> filter(group==2) |> select(prtvede2) |> table() / 601
countrydf |> filter(group==3) |> select(prtvede2) |> table() / 562

262+232+4+3+4+57

countrydf |> select(prtvede2) |> table() / 2203

countrydf |> filter(group==3)  |> nrow()
countrydf |> select(group) |> table()

coherence <- function(v = rep(0,5), R = diag(rep(1,length(v)))) c(0.5*v%*%(as.matrix(R)-diag(rep(1,length(v))))%*%v)
Rs <- R |> select(-item) |> nest(.by = "group") |> 
  mutate(vectors = data |> map(function(x) eigen(x)$vectors),
         values = data |> map(function(x) eigen(x)$values),
         val1 = values |> map_dbl(function(x) x[1]),
         ev1 = vectors |> map(function(x) x[,1] * sign(x[3,1])),
         ev1_face = ev1 |> map(function(v) v/max(abs(v))),
         ev1_corner = ev1 |> map(function(v) sign(v)),
         val2 = values |> map_dbl(function(x) x[2]),
         ev2 = vectors |> map(function(x) x[,2] * sign(x[3,2])),
         coh_face = map2_dbl(ev1_face, data, coherence),
         coh_corner = map2_dbl(ev1_corner, data, coherence),
         percent_variance = val1/5)
Rs

