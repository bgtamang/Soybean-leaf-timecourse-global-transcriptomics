# Deploy to shinyapps.io
library(rsconnect)

# Set your shinyapps.io credentials before deploying:
# rsconnect::setAccountInfo(name = 'bgtamang', token = 'YOUR_TOKEN', secret = 'YOUR_SECRET')
# Or run rsconnect::setAccountInfo(...) in your R console first.

cat("Account configured. Deploying...\n")

rsconnect::deployApp(
  appDir = "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase3-Refined-Analysis/GmJAG1_Dashboard",
  appName = "GmJAG1-Dashboard",
  account = "bgtamang",
  forceUpdate = TRUE
)
