library(IsoCorrectoR)
IsoCorrection(MeasurementFile = "./MeasurementFile.csv", ElementFile = "./ElementFile.csv", MoleculeFile = "./MoleculeFile.csv",
       CorrectTracerImpurity = FALSE, CorrectTracerElementCore = TRUE,
       CalculateMeanEnrichment = TRUE, UltraHighRes = FALSE, DirOut = ".",
       FileOut = "result", FileOutFormat = "csv", ReturnResultsObject = TRUE,
       CorrectAlsoMonoisotopic = FALSE, CalculationThreshold = 10^-8,
       CalculationThreshold_UHR = 8, verbose = FALSE, Testmode = FALSE)
