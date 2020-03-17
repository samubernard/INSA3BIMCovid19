# INSA3BIMCovid19
Spatial epidemiological model for covid-19 epidemics in France at INSA Lyon

- SEIR epidemiological model
- Spatial diffusion for local spread
- Indivual-based-model layer for travellers and 'super-infectors'

Effect of
- No isolation measure
- mild isolation measures: restricting travel, avoiding crowds
- strong measures: stay at home except for work and grocery shopping
- emergency measures: stay at home, no outside work except emergency services

How and when to lift these measures?

## About Covid-19

From https://www.who.int 

> Coronaviruses are a large family of viruses which may cause illness in animals or humans.  In humans, several coronaviruses are known to cause respiratory infections ranging from the common cold to more severe diseases such as Middle East Respiratory Syndrome (MERS) and Severe Acute Respiratory Syndrome (SARS). The most recently discovered coronavirus causes coronavirus disease COVID-19.

> COVID-19 is the infectious disease caused by the most recently discovered coronavirus. This new virus and disease were unknown before the outbreak began in Wuhan, China, in December 2019.

> The most common symptoms of COVID-19 are fever, tiredness, and dry cough. Some patients may have aches and pains, nasal congestion, runny nose, sore throat or diarrhea. These symptoms are usually mild and begin gradually. Some people become infected but don’t develop any symptoms and don't feel unwell. Most people (about 80%) recover from the disease without needing special treatment. Around 1 out of every 6 people who gets COVID-19 becomes seriously ill and develops difficulty breathing. Older people, and those with underlying medical problems like high blood pressure, heart problems or diabetes, are more likely to develop serious illness. People with fever, cough and difficulty breathing should seek medical attention.

> People can catch COVID-19 from others who have the virus. The disease can spread from person to person through small droplets from the nose or mouth which are spread when a person with COVID-19 coughs or exhales. These droplets land on objects and surfaces around the person. Other people then catch COVID-19 by touching these objects or surfaces, then touching their eyes, nose or mouth. People can also catch COVID-19 if they breathe in droplets from a person with COVID-19 who coughs out or exhales droplets. This is why it is important to stay more than 1 meter (3 feet) away from a person who is sick.
WHO is assessing ongoing research on the ways COVID-19 is spread and will continue to share updated findings.    

> Studies to date suggest that the virus that causes COVID-19 is mainly transmitted through contact with respiratory droplets rather than through the air.  See previous answer on “How does COVID-19 spread?”

> The main way the disease spreads is through respiratory droplets expelled by someone who is coughing. The risk of catching COVID-19 from someone with no symptoms at all is very low. However, many people with COVID-19 experience only mild symptoms. This is particularly true at the early stages of the disease. It is therefore possible to catch COVID-19 from someone who has, for example, just a mild cough and does not feel ill.  WHO is assessing ongoing research on the period of transmission of COVID-19 and will continue to share updated findings.    


### An excellent overview of the main issues we are facing:

+ https://medium.com/@tomaspueyo/coronavirus-act-today-or-people-will-die-f4d3d9cd99ca

+ The Korean Cluster https://graphics.reuters.com/CHINA-HEALTH-SOUTHKOREA-CLUSTERS/0100B5G33SB/index.html

+ How much is coronavirus spreading under the radar? https://www.nature.com/articles/d41586-020-00760-8]]Datasets
https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases

### R0 reproductive number estimates
We estimated that the median daily reproduction number (Rt) in Wuhan declined from 2·35 (95% CI 1·15–4·77) 1 week before travel restrictions were introduced on Jan 23, 2020, to 1·05 (0·41–2·39) 1 week after. Kucharski et al 2020 The Lancet

*Disease progression times* 	

Type                          time
----------------              --------
Incubation period 	          5.2 days
infectious period 	          2.9 days
delay onset-to-confirmation 	6.1 days


SEIR Model
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model


