

do_everything("2020-2-3")


#####

test <- get_all() %>%
  check_all() %>%
  join_data()

test2 <- test %>%
  dplyr::filter(lubridate::year(start) == 2019)


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
