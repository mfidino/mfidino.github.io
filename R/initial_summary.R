

do_everything("2021-03-01")

# For year end ARCs reporting
arcs <- get_all() %>%
  check_all() %>%
  join_data()


str(dat)

dat$start <-lubridate::ymd(dat$start)
dat$end <-lubridate::ymd(dat$end)

# add more time
test <- dplyr::bind_rows(
  dat,
  data.frame(
    start = seq(tail(dat$start,1), length.out = 50, by = "14 days"),
    end = seq(tail(dat$end,1), length.out = 50, by = "14 days"),
    stringsAsFactors = FALSE
  )
)

dat <- apply(dat, 2, function(x) lubridate::ymd(as.character(x)))
#####

test <- get_all() %>%
  check_all() %>%
  join_data()
test2 <- test %>%
  dplyr::filter(lubridate::year(start) == 2020)

write.csv(test, "./data/meeting_days.csv", row.names = FALSE)


  dplyr::filter(test$end < lubridate::ymd_hms('2020-1-1 0:00:00'))

first_pass <- test2 %>%
  dplyr::mutate(dur = as.numeric(.$end - .$start)) %>%
  dplyr::mutate(prop_time = dur / sum(dur))

x <- first_pass[,c("project_name", "client_category")] %>%
  dplyr::distinct()

by_project <- first_pass %>%
  dplyr::group_by(project_name) %>%
  dplyr::summarise(prop_time = sum(prop_time)) %>%
  dplyr::arrange(dplyr::desc(prop_time)) %>%
  dplyr::left_join(., x, by = "project_name") %>%
  dplyr::select(project_name, client_category, prop_time) %>%
  dplyr::distinct() %>%
  na.omit()

by_project$prop_time <- by_project$prop_time / sum(by_project$prop_time)

by_project$prop_time <- by_project$prop_time * 100

by_project$prop_time <- round(by_project$prop_time, 2)
yo <- data.frame(by_project[order(by_project$prop_time, decreasing = TRUE),])

yo$used<- 0
yo$used[7] <- 1

# field work stuff
sum(yo$prop_time[c(48,13)])
yo$used[c(13,48)] <- 1

# bats
yo$used[37] <- 1

# uwin
sum(yo$prop_time[yo$client_category == "UWIN"]) + 6.49 + 24.30
yo$used[yo$client_category == "UWIN"] <- 1
yo$used[1] <- 1
yo$used[4] <- 1



yo$used[3] <- 1



sum(yo$prop_time[grep("rat", yo$project_name)])
yo$used[grep("rat", yo$project_name)] <- 1

yo$used[28] <- 1

yo$used[34] <- 1


51.67 + 1.04


sum(yo$prop_time[c(11, 20,21,23,27,31,41)])

yo$used[c(11, 20,21,23,27,31,41)] <- 1

yo[!yo$used,]
