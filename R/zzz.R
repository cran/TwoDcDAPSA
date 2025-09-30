.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      "Thank you for using TwoDcDAPSA (PJC Calculator).",
      "Ranges used for cDAPSA inputs:",
      " - Swollen Joint Count (0-66)",
      " - Tender Joint Count  (0-68)",
      " - Pain                (0-10)",
      " - Patient Global      (0-10)",
      sep = "\n"
    )
  )
}
