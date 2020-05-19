

do_everything("2020-5-11")


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
  dplyr::filter(lubridate::year(start) == 2019)

write.csv(test, "./data/meeting_days.csv", row.names = FALSE)


  dplyr::filter(test$end < lubridate::ymd_hms('2020-1-1 0:00:00'))

first_pass <- test2 %>%
  dplyr::mutate(dur = as.numeric(.$end - .$start)) %>%
  dplyr::mutate(prop_time = dur / sum(dur))

by_project <- first_pass %>%
  dplyr::group_by(project_name) %>%
  dplyr::summarise(prop_time = sum(prop_time)) %>%
  dplyr::arrange(dplyr::desc(prop_time)) %>%
  dplyr::left_join(., x, by = "project_name") %>%
  dplyr::select(project_name, client_category, prop_time) %>%
  dplyr::distinct() %>%
  na.omit()
