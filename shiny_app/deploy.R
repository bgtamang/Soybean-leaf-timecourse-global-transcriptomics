# Deploy to shinyapps.io
library(rsconnect)

# Set your shinyapps.io credentials before deploying:
# rsconnect::setAccountInfo(name = 'bgtamang', token = 'YOUR_TOKEN', secret = 'YOUR_SECRET')
# Or run rsconnect::setAccountInfo(...) in your R console first.

cat("Account configured. Deploying...\n")

rsconnect::deployApp(
  appDir = ".",
  appName = "GmJAG1-Dashboard",
  account = "bgtamang",
  forceUpdate = TRUE
)
