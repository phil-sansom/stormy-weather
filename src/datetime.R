extract.dates = function(t0, times, calendar, units) {
  
  dd = addtimes(t0, times, calendar, units)
  dd = substr(dd, 1, 11)
  return(dd)
  
}

difftimes = function(date1, date2, calendar = "standard") {
  
  cal = tolower(calendar)
  if (cal %in% c("gregorian","standard","proleptic_gregorian")) {
    d1 = as.POSIXlt(date1, format = "%Y%m%dT%H%M%S", tz = "UTC")
    d2 = as.POSIXlt(date2, format = "%Y%m%dT%H%M%S", tz = "UTC")
    diff = as.numeric(difftime(d1, d2, units = "hours"))
  } else if (cal %in% c("noleap","365_day","all_leap","366_day","360_day")) {
    d1 = to.POSIXct(date1, calendar)
    d2 = to.POSIXct(date2, calendar)
    diff = (d1 - d2)/60/60
  } else {
    diff = NA
  }
  return(diff)
  
}

addtimes = function(date, diff, calendar = "standard", units = "hours") {
  
  cal = tolower(calendar)
  if (cal %in% c("gregorian","standard","proleptic_gregorian")) {
    t0 = as.POSIXlt(date, format = "%Y%m%dT%H%M%S")
    td = as.difftime(diff, units = units)
    dd = format(t0 + td, "%Y%m%dT%H%M%S", tz = "UTC")
  } else if (cal %in% c("noleap","365_day","all_leap","366_day","360_day")) {
    t0 = to.POSIXct(date, calendar)
    td = to.difftime(diff, units)
    dd = to.character(t0 + td, calendar)
  } else {
    dd = NA
  }
  return(dd)
  
}



to.POSIXct = function(x, calendar) {
  
  ## Seconds since 1970-01-01 00:00:00
  yyyy = as.numeric(substr(x, 1, 4))
  mm   = as.numeric(substr(x, 5, 6))
  dd   = as.numeric(substr(x, 7, 8))
  HH   = as.numeric(substr(x, 10, 11))
  MM   = as.numeric(substr(x, 12, 13))
  SS   = as.numeric(substr(x, 14, 15))
  if (is.na(HH))
    HH = 0
  if (is.na(MM))
    MM = 0
  if (is.na(SS))
    SS = 0
  cal = tolower(calendar)
  if (cal %in% c("noleap","365_day")) {
    mlen = c(0,31,59,90,120,151,181,212,243,273,304,334)
    y = (yyyy - 1970)*365*24*60*60 + mlen[mm]*24*60*60 + 
      (dd - 1)*24*60*60 + HH*60*60 + MM*60 + SS
  } else if (cal %in% c("all_leap","366_day")){
    mlen = c(0,31,60,91,121,152,182,213,244,274,305,335)
    y = (yyyy - 1970)*366*24*60*60 + mlen[mm]*24*60*60 + 
      (dd - 1)*24*60*60 + HH*60*60 + MM*60 + SS
  } else if (cal == "360_day") {
    y = (yyyy - 1970)*12*30*24*60*60 + (mm - 1)*30*24*60*60 + 
      (dd - 1)*24*60*60 + HH*60*60 + MM*60 + SS
  } else {
    y = NA
  }
  return(y)
  
}

to.difftime = function(x, units) {
  
  u = tolower(units)
  if (u %in% c("days","day","ds","d")) {
    y = x*24*60*60
  } else if (u %in% c("hours","hour","hrs","hr","hs","h")) {
    y = x*60*60
  } else if (u %in% c("minutes","minute","mins","min")) {
    y = x*60
  } else if (u %in% c("secs","sec","s")) {
    y = x
  } else {
    y = NA
  }
  
  return(y)
  
}

to.character = function(x, calendar) {
  
  cal = tolower(calendar)
  if (cal %in% c("noleap","365_day")) {
    yyyy = 1970 + x %/% (365*24*60*60)
    xx = x %% (365*24*60*60)
    mlen = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
    m = xx %/% (24*60*60)
    mm = numeric(length(x))
    dd = numeric(length(x))
    for (i in 1:12) {
      mask = mlen[i] <= m & m < mlen[i+1]
      mm[mask] = i
      dd[mask] = m[mask] - mlen[i] + 1
    }
  } else if (cal %in% c("all_leap","366_day")){
    yyyy = 1970 + x %/% (366*24*60*60)
    xx = x %% (366*24*60*60)
    mlen = c(0,31,60,91,121,152,182,213,244,274,305,335,366)
    m = xx %/% (24*60*60)
    mm = numeric(length(x))
    dd = numeric(length(x))
    for (i in 1:12) {
      mask = mlen[i] <= m & m < mlen[i+1]
      mm[mask] = i
      dd[mask] = m[mask] - mlen[i] + 1
    }
  } else if (cal == "360_day") {
    yyyy = 1970 + x %/% (360*24*60*60)
    xx = x  %%  (360*24*60*60)
    mm = xx %/% (30*24*60*60) + 1
    xx = xx %%  (30*24*60*60)
    dd = xx %/% (24*60*60) + 1
  } else {
    xx = NA
  }
  xx = xx %%  (24*60*60)
  HH = xx %/% (60*60)
  xx = xx %%  (60*60)
  MM = xx %/% 60
  xx = xx %%  60
  SS = xx
  
  yyyy = sprintf("%04i", yyyy)
  mm = sprintf("%02i", mm)
  dd = sprintf("%02i", dd)
  HH = sprintf("%02i", HH)
  MM = sprintf("%02i", MM)
  SS = sprintf("%02i", round(SS))
  
  z = paste0(yyyy,mm,dd,"T",HH,MM,SS)

  return(z)
  
}
