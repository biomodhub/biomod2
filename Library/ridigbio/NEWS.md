ridigbio 0.4.1
===========
### Fixes
* Styling Fixes and Meta Fields fixes and tests 

ridigbio 0.4.0
===========
### Enhancements
* Offline build capability added

ridigbio 0.3.9
===========

### Documentation
* Pkgdown site created!

ridigbio 0.3.8
===========

### NEW FEATURES
* Default fields returned to users have been updated to return research-grade fields. Previously, we returned the datecollected field by default, which we do not recommend to be used in scientific research.
* What is datecollected? This is a field created during the data ingestion process. When a data provider does not provide a full date in the Darwin Core eventDate field, this complete value or the missing parts (i.e., month and/or day) are randomly generated and thus may lack any real meaning. The generated dates are difficult to detect, as they are randomly distributed. We are currently working to modify our ingestion pipeline to avoid randomly generating dates. However, dates remain an issue across biodiversity aggregators and the solution is not clear (see GBIF for example).
* To prevent user misuse of this term, we will no longer be providing the datecollected field by default and will instead be returning the following fields: c("data.dwc:eventDate", "data.dwc:year",  "data.dwc:month", "data.dwc:day").
* Please be advised that these fields are not in date format. Instead, all dates will be text strings. There are many ways to convert these to dates, for example, see gatoRs remove_duplicate function or ridigbio proposed solution here.

### Documentation
* Additional function documentation has been added to streamline data use.
