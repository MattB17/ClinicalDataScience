library(tidyverse)
library(magrittr)
library(bigrquery)
library(caret)

conn <- DBI::dbConnect(drv=bigquery(),
                       project="learnclinicaldatascience")

# getStats is a course provided function
getStats <- function(df, ...) {
  df %>%
    select_(.dots = lazyeval::lazy_dots(...)) %>%
    mutate_all(funs(factor(., levels=c(1,0)))) %>%
    table() %>%
    confusionMatrix()
}

hypertension <- tbl(conn, "course3_data.hypertension_goldstandard")
diagnoses_icd <- tbl(conn, "mimic3_demo.DIAGNOSES_ICD")

# returns diagnoses based on icd codes
getIcdData <- function(icdCodes) {
  diagnoses_icd %>%
    filter(ICD9_CODE %in% icdCodes) %>%
    distinct(SUBJECT_ID) %>%
    mutate(predictor = 1)
}

# joins joinDf to hypertension in order to use it as a predictor variable
# then analyzes the result
performJoinAndGetStats <- function(joinDf) {
  hypertension %>%
    left_join(joinDf) %>%
    mutate(predictor = coalesce(predictor, 0)) %>%
    collect() %>%
    getStats(predictor, HYPERTENSION)
}

## joins the two dataframes and creates an indicator variable for
## each patient that experiences an event above threshVal
getThreshFromInnerJoin <- function(mainDf, newDf, threshVal) {
  mainDf %>%
    inner_join(newDf) %>%
    mutate(counter=1) %>%
    group_by(SUBJECT_ID) %>%
    summarize(countVar = sum(counter, na.rm=TRUE)) %>%
    mutate(predictor = case_when(countVar > threshVal ~ 1, TRUE ~ 0))
}

## to get number of unique subjects after inner join
performInnerJoinCountDistinct <- function(mainDf, joinDf) {
  mainDf %>%
    inner_join(joinDf) %>%
    distinct(SUBJECT_ID) %>%
    mutate(predictor = 1)
}

# creates an indicator for each patient if the proportion of prescriptions
# for hypertensive are at or above threshVal
getHypertensivePrescriptionsAboveThresh <- function(threshVal) {
  antihypertensives <- tbl(conn, "course3_data.D_ANTIHYPERTENSIVES") %>%
    mutate(hypertension_prescription_counter = 1)
  
  prescriptions %>%
    left_join(antihypertensives) %>%
    mutate(hypertension_prescription_counter = 
             coalesce(hypertension_prescription_counter, 0),
           prescription_counter = 1) %>%
    group_by(SUBJECT_ID) %>%
    summarize(hypertension_prescription_count = 
                sum(hypertension_prescription_counter, na.rm=TRUE),
              prescription_count = 
                sum(prescription_counter, na.rm=TRUE)) %>%
    mutate(predictor = 
             case_when((hypertension_prescription_count / prescription_count) >= threshVal ~ 1,
                       TRUE ~ 0))
}

## combines two indicators into one and gets statistics on predictive
## power of the resulting indicator
getStatsForTwoIndicators <- function(df1, df2) {
  hypertension %>%
    left_join(df1) %>%
    left_join(df2) %>%
    mutate(predictor2 = coalesce(predictor2, 0),
           predictor = coalesce(predictor, 0)) %>%
    mutate(ind_var = case_when((predictor2 == 1) && (predictor == 1) ~ 1,
                               TRUE ~ 0)) %>%
    collect() %>%
    getStats(ind_var, HYPERTENSION)
}

# ICD 401.0 only
icdData <- getIcdData(c("4010"))
performJoinAndGetStats(icdData)

# ICD 572.3 only
icdData <- getIcdData(c("5723"))
performJoinAndGetStats(icdData) 

# Atenolol prescription only
prescriptions <- tbl(conn, "mimic3_demo.PRESCRIPTIONS")

atenolol <- prescriptions %>%
  filter(tolower(DRUG) == 'atenolol') %>%
  distinct(SUBJECT_ID) %>%
  mutate(predictor = 1)

performJoinAndGetStats(atenolol)

# At least two essential hypertension codes
essential_hypertension_thresh <- diagnoses_icd %>%
  mutate(essential_hypertension_counter = 
           case_when(ICD9_CODE %in% c("4010", "4011", "4019") ~ 1,
                     TRUE ~ 0)) %>%
  group_by(SUBJECT_ID) %>%
  summarize(essential_hypertension_count = 
              sum(essential_hypertension_counter, na.rm=TRUE)) %>%
  mutate(predictor = 
           case_when(essential_hypertension_count > 1 ~ 1,
                     TRUE ~ 0))

performJoinAndGetStats(essential_hypertension_thresh)

# At least two non essential hypertension codes
non_essential_hypertension_codes <- tbl(conn, "mimic3_demo.D_ICD_DIAGNOSES") %>%
  filter(tolower(LONG_TITLE) %like% '%hypertension%',
         !ICD9_CODE %in% c("4010", "4011", "4019")) %>%
  select(ICD9_CODE)

non_essential_hypertension_thresh <- getThreshFromInnerJoin(
  diagnoses_icd, non_essential_hypertension_codes)

performJoinAndGetStats(non_essential_hypertension_thresh)

# At least 2 hypertension prescriptions
antihypertensives <- tbl(conn, "course3_data.D_ANTIHYPERTENSIVES")

hypertension_prescriptions_thresh <- getThreshFromInnerJoin(
  prescriptions, antihypertensives, 1)

performJoinAndGetStats(hypertension_prescriptions_thresh)

# At least 10% of prescriptions are for hypertension
hypertension_prescriptions_thresh <- getHypertensivePrescriptionsAboveThresh(0.1)

performJoinAndGetStats(hypertension_prescriptions_thresh)

# essential hypertension icd
essential_hypertension_icd <- getIcdData(c("4010", "4011", "4019"))

performJoinAndGetStats(essential_hypertension_icd)

### any hypertensive icd
hypertensive_codes <- tbl(conn, "mimic3_demo.D_ICD_DIAGNOSES") %>%
  select(ICD9_CODE)

hypertension_icd <- performInnerJoinCountDistinct(
  diagnoses_icd, hypertensive_codes)

performJoinAndGetStats(hypertension_icd)

### any anti hypertensive prescription
antihypertensives <- tbl(conn, "course3_data.D_ANTIHYPERTENSIVES")

anti_hypertensive_prescription <- performInnerJoinCountDistinct(
  prescriptions, antihypertensives)

performJoinAndGetStats(anti_hypertensive_prescription)

### essential hypertension icd and any anti hypertensive prescription
anti_hypertensive_prescription <- performInnerJoinCountDistinct(
  prescriptions, antihypertensives)

essential_hypertension_icd <- getIcdData(c("4010", "4011", "4019")) %>%
  mutate(predictor2 = predictor) %>%
  select(SUBJECT_ID, predictor2)

getStatsForTwoIndicators(anti_hypertensive_prescription, essential_hypertension_icd)

### essential hypertension icd and at least 10% anti hypertensive prescription
hypertension_prescriptions_thresh <- getHypertensivePrescriptionsAboveThresh(0.1)

getStatsForTwoIndicators(hypertension_prescriptions_thresh, essential_hypertension_icd)