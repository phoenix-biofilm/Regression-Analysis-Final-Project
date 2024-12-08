## Final Project Data Clearning

## Packages
if (!require("tidyverse")){
  install.packages("tidyverse")
}
library(tidyverse)

## Open and view raw data
gun_violence <- read.csv("C:/Users/willphoe/Documents/Regression Analysis Final Project/Gun Violence Data.csv", header = TRUE)

head(gun_violence)
str(gun_violence)
summary(gun_violence)

## Remove unnecessary variables & rename variables appropriately
gun_violence_clean <- subset(gun_violence, select = c(incident_id, date, state, 
                                                city_or_county, n_killed, n_injured, 
                                                congressional_district, gun_stolen, 
                                                gun_type, latitude, longitude, 
                                                n_guns_involved, participant_age, 
                                                participant_gender, participant_status, 
                                                participant_type, state_house_district, 
                                                state_senate_district))

## Remove observations with nan values
gun_violence_clean <- gun_violence_clean |> na.omit()

## Change variable types as necessary
gun_violence_clean <-
  gun_violence_clean |> transform(state = factor(state),
                                  city_or_county = factor(city_or_county))


## Reformat gun_stolen, gun_type, participant_age, participant_gender, participant_status, 
## and participant_types
nobs <- nrow(gun_violence_clean)

# gun_stolen
gun_violence_clean$gun_stolen <- strsplit(gun_violence_clean$gun_stolen,"[||]|[|]")

# gun_type
gun_violence_clean$gun_type <- strsplit(gun_violence_clean$gun_type,"[||]|[|]")

# participant_age
gun_violence_clean$participant_age <- strsplit(gun_violence_clean$participant_age,"[||]|[|]")

# participant_gender
gun_violence_clean$participant_gender <- strsplit(gun_violence_clean$participant_gender,
                                                               "[||]|[|]")

# participant_status
gun_violence_clean$participant_status <- strsplit(gun_violence_clean$participant_status,
                                                               "[||]|[|]")

# participant_type
gun_violence_clean$participant_type <- strsplit(gun_violence_clean$participant_type,
                                                             "[||]|[|]")

## Create aggregate columns
gun_violence_clean$n_guns_stolen <- integer(length=nobs)
gun_violence_clean$n_guns_not_stolen <- integer(length=nobs)
gun_violence_clean$n_guns_unknown <- integer(length=nobs) 
gun_violence_clean$n_guns_type_unknown <- integer(length=nobs)
gun_violence_clean$n_handguns <- integer(length=nobs)
gun_violence_clean$n_shotguns <- integer(length=nobs)
gun_violence_clean$n_gun_type_known <- integer(length=nobs)
gun_violence_clean$n_agel20 <- integer(length=nobs)
gun_violence_clean$n_age20_30 <- integer(length=nobs)
gun_violence_clean$n_age30_40 <- integer(length=nobs)
gun_violence_clean$n_age40_50 <- integer(length=nobs)
gun_violence_clean$n_age50_60 <- integer(length=nobs)
gun_violence_clean$n_age60_70 <- integer(length=nobs)
gun_violence_clean$n_ageg70 <- integer(length=nobs)
gun_violence_clean$n_male <- integer(length=nobs)
gun_violence_clean$n_female <- integer(length=nobs)
gun_violence_clean$n_arrested <- integer(length=nobs)
gun_violence_clean$n_unharmed <- integer(length=nobs)
gun_violence_clean$n_suspects <- integer(length=nobs)
gun_violence_clean$n_victims <- integer(length=nobs)
gun_violence_clean$n_killed_div_n_injured <- numeric(length=nobs)


## Remove blank entries in lists
for (i in 1:nobs){
  gun_violence_clean$gun_stolen[[i]] <- Filter(function(x)x!="", gun_violence_clean$gun_stolen[[i]])
  gun_violence_clean$gun_type[[i]] <- Filter(function(x)x!="", gun_violence_clean$gun_type[[i]])
  gun_violence_clean$participant_age[[i]] <- Filter(function(x)x!="", gun_violence_clean$participant_age[[i]])
  gun_violence_clean$participant_gender[[i]] <- Filter(function(x)x!="", gun_violence_clean$participant_gender[[i]])
  gun_violence_clean$participant_status[[i]] <- Filter(function(x)x!="", gun_violence_clean$participant_status[[i]])
  gun_violence_clean$participant_type[[i]] <- Filter(function(x)x!="", gun_violence_clean$participant_type[[i]])
}
print("done 1")

## Remove ID numbers & calculate aggregate counts
for (i in 1:nobs){
  for (j in 1:length(gun_violence_clean$gun_stolen[[i]])){
    gun_violence_clean$gun_stolen[[i]][j] <- gsub(".*::", "", gun_violence_clean$gun_stolen[[i]][j])
    gun_violence_clean$gun_stolen[[i]][j] <- gsub(".*:", "", gun_violence_clean$gun_stolen[[i]][j])
    if (is.na(gun_violence_clean$gun_stolen[[i]][j])){
      break
    } else if (gun_violence_clean$gun_stolen[[i]][j] == "Stolen"){
      gun_violence_clean$n_guns_stolen[i] <- gun_violence_clean$n_guns_stolen[i] + 1
    } else if (gun_violence_clean$gun_stolen[[i]][j] == "Not-stolen"){
      gun_violence_clean$n_guns_not_stolen[i] <- gun_violence_clean$n_guns_not_stolen[i] + 1
    } else {
      gun_violence_clean$n_guns_unknown[i] <- gun_violence_clean$n_guns_unknown[i] + 1
    }
  }
  for (j in 1:length(gun_violence_clean$gun_type[[i]])){
    gun_violence_clean$gun_type[[i]][j] <- gsub(".*::", "", gun_violence_clean$gun_type[[i]][j])
    gun_violence_clean$gun_type[[i]][j] <- gsub(".*:", "", gun_violence_clean$gun_type[[i]][j])
    if (is.na(gun_violence_clean$gun_type[[i]][j])){
      break
    } else if(gun_violence_clean$gun_type[[i]][j] == "Unknown"){
      gun_violence_clean$n_guns_type_unknown[i] <- gun_violence_clean$n_guns_type_unknown[i]+ 1
    } else if (gun_violence_clean$gun_type[[i]][j] == "Handgun"){
      gun_violence_clean$n_handguns[i] <- gun_violence_clean$n_handguns[i] + 1
    } else if (gun_violence_clean$gun_type[[i]][j] == "Shotgun"){
      gun_violence_clean$n_shotguns[i] <- gun_violence_clean$n_shotguns[i] + 1
    } else {
      gun_violence_clean$n_guns_type_known[i] <- gun_violence_clean$n_guns_type_known[i] + 1
    }
  }
  for (j in 1:length(gun_violence_clean$participant_age[[i]])){
    gun_violence_clean$participant_age[[i]][j] <- gsub(".*::", "", gun_violence_clean$participant_age[[i]][j])
    gun_violence_clean$participant_age[[i]][j] <- gsub(".*:", "", gun_violence_clean$participant_age[[i]][j])
    if (is.na(gun_violence_clean$participant_age[[i]][j])){
      break
    } else if (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 20){
      gun_violence_clean$n_agel20[i] <- gun_violence_clean$n_agel20[i] + 1
    } else if ((20 <= as.numeric(gun_violence_clean$participant_age[[i]][j])) & 
               (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 30)){
      gun_violence_clean$n_age20_30[i] <- gun_violence_clean$n_age20_30[i] + 1
    } else if ((30 <= as.numeric(gun_violence_clean$participant_age[[i]][j])) & 
               (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 40)){
      gun_violence_clean$n_age30_40[i] <- gun_violence_clean$n_age30_40[i] + 1
    } else if ((40 <= as.numeric(gun_violence_clean$participant_age[[i]][j])) & 
               (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 50)){
      gun_violence_clean$n_age40_50[i] <- gun_violence_clean$n_age40_50[i] + 1
    } else if ((50 <= as.numeric(gun_violence_clean$participant_age[[i]][j])) & 
               (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 60)){
      gun_violence_clean$n_age50_60[i] <- gun_violence_clean$n_age50_60[i] + 1
    } else if ((60 <= as.numeric(gun_violence_clean$participant_age[[i]][j])) & 
               (as.numeric(gun_violence_clean$participant_age[[i]][j]) < 70)){
      gun_violence_clean$n_age60_70[i] <- gun_violence_clean$n_age60_70[i] + 1
    } else if (70 < as.numeric(gun_violence_clean$participant_age[[i]][j])){
      gun_violence_clean$n_ageg70[i] <- gun_violence_clean$n_ageg70[i] + 1
    } else {
      break
    }
  }
  for (j in 1:length(gun_violence_clean$participant_gender[[i]])){
    gun_violence_clean$participant_gender[[i]][j] <- gsub(".*::", "", gun_violence_clean$participant_gender[[i]][j])
    gun_violence_clean$participant_gender[[i]][j] <- gsub(".*:", "", gun_violence_clean$participant_gender[[i]][j])
    if (is.na(gun_violence_clean$participant_gender[[i]][j])){
      break
    } else if (gun_violence_clean$participant_gender[[i]][j] == "Female"){
      gun_violence_clean$n_female[i] <- gun_violence_clean$n_female[i] + 1
    } else if (gun_violence_clean$participant_gender[[i]][j] == "Male"){
      gun_violence_clean$n_male[i] <- gun_violence_clean$n_male[i] + 1
    } else {
      break
    }
  }
  for (j in 1:length(gun_violence_clean$participant_status[[i]])){
    gun_violence_clean$participant_status[[i]][j] <- gsub(".*::", "", gun_violence_clean$participant_status[[i]][j])
    gun_violence_clean$participant_status[[i]][j] <- gsub(".*:", "", gun_violence_clean$participant_status[[i]][j])
    if (is.na(gun_violence_clean$participant_status[[i]][j])){
      break
    } else if (gun_violence_clean$participant_status[[i]][j] == "Arrested"){
      gun_violence_clean$n_arrested[i] <- gun_violence_clean$n_arrested[i] + 1
    } else if (gun_violence_clean$participant_status[[i]][j] == "Unharmed"){
      gun_violence_clean$n_unharmed[i] <- gun_violence_clean$n_unharmed[i] + 1
    } else {
      break
    }
  }
  for (j in 1:length(gun_violence_clean$participant_type[[i]])){
    gun_violence_clean$participant_type[[i]][j] <- gsub(".*::", "", gun_violence_clean$participant_type[[i]][j])
    gun_violence_clean$participant_type[[i]][j] <- gsub(".*:", "", gun_violence_clean$participant_type[[i]][j])
    if (is.na(gun_violence_clean$participant_type[[i]][j])){
      break
    } else if (gun_violence_clean$participant_type[[i]][j] == "Subject-Suspect"){
      gun_violence_clean$n_suspects[i] <- gun_violence_clean$n_suspects[i] + 1
    } else if (gun_violence_clean$participant_type[[i]][j] == "Victim"){
      gun_violence_clean$n_victims[i] <- gun_violence_clean$n_victims[i] + 1
    } else {
      break
    }
  }
  
}
print("Done 1")
gun_violence_clean$n_killed_div_n_injured <- gun_violence_clean$n_killed/gun_violence_clean$n_injured


## Delete original variables
gun_violence_clean <- subset(gun_violence_clean, select = c(incident_id, date, 
                                                            state, city_or_county,
                                                            n_killed, n_injured, 
                                                            n_killed_div_n_injured, 
                                                            congressional_district,
                                                            n_guns_stolen, n_guns_not_stolen,
                                                            n_guns_unknown, n_guns_type_unknown,
                                                            n_handguns, n_shotguns,
                                                            n_guns_type_known, latitude,
                                                            longitude, n_guns_involved,
                                                            n_agel20, n_age20_30, 
                                                            n_age30_40, n_age40_50, 
                                                            n_age50_60, n_age60_70, 
                                                            n_ageg70, n_male, n_female,
                                                            n_arrested, n_unharmed, 
                                                            n_suspects, n_victims, 
                                                            state_house_district, 
                                                            state_senate_district))

print("done 4")

## Remove observations with nan values again
gun_violence_clean <- gun_violence_clean |> na.omit()

## Review cleaned data, including summaries and visual representations for each predictor variable
head(gun_violence_clean)
str(gun_violence_clean)
summary(gun_violence_clean)



## Save clean data to a csv file
write.csv(gun_violence_clean, "gun_violence_clean.csv", row.names = TRUE)

