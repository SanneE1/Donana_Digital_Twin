

dry_wet_years <- function (r_temp, r_precip) {
  
  a <- global(r_temp, mean, na.rm = T)
  a$year <- year(terra::time(r_temp))
  a$month <- month(terra::time(r_temp))
  a <- a %>%
    group_by(month) %>%
    mutate(z_temp = scale(mean),
           .keep = "unused")
  
  b <- global(r_precip, mean, na.rm = T)
  b$year <- year(terra::time(r_precip))
  b$month <- month(terra::time(r_precip))
  b <- b %>%
    group_by(month) %>%
    mutate(z_precip = scale(mean),
           .keep = "unused")
  
  c <- left_join(a,b)
 
  return(c) 
}