# --- 000_layout ---------------------------------------------------------------

# --- Figures 1 and S1 ---------------------------------------------------------

fs_images <- tribble(
  ~fig,    ~panel,        ~x,     ~y,     ~w,     ~h,
  "fig1",  "Alluses",     210.45,   0.00, 172.94, 172.95,
  "fig1",  "Fisheries",     0.00, 178.77, 170.37, 170.39,
  "fig1",  "Aquarium",    210.45, 178.77, 170.37, 170.39,
  "fig1",  "Aquaculture",   0.00, 348.61, 170.37, 170.39,
  "fig1",  "Gamefish",    210.45, 348.61, 170.37, 170.39,
  "figS1", "Alluses",     210.11,   0.00, 173.75, 174.05,
  "figS1", "Fisheries",     0.00, 178.77, 170.10, 170.39,
  "figS1", "Aquarium",    210.11, 178.77, 170.10, 170.39,
  "figS1", "Aquaculture",   0.00, 348.61, 170.10, 170.39,
  "figS1", "Gamefish",    210.11, 348.61, 170.10, 170.39
)

circle_place <- list(
  fig1  = list(box = c(-0.70, 7.99, 172.68, 168.94), clip = c(3.25, 13.95, 170.37, 172.38)),
  figS1 = list(box = c(1.50, 11.95, 168.36, 168.61), clip = NULL)
)

fig1_masks <- tribble(
  ~x, ~y,
  list(c(227.58, 233.59, 233.59, 227.58)), list(c(349.69, 349.69, 354.64, 354.64)),
  list(c(229.58, 231.04, 231.04, 229.58)), list(c(364.46, 364.46, 367.12, 367.12)),
  list(c(288.87, 309.26, 309.26, 288.87)), list(c(508.79, 508.79, 512.42, 512.42)),
  list(c(212.95, 218.93, 218.93, 212.95)), list(c(421.63, 421.63, 456.42, 456.42)),
  list(c(78.45, 100.09, 100.09, 78.45)), list(c(337.65, 337.65, 346.37, 346.37)),
  list(c(1.93, 7.97, 7.97, 1.93)), list(c(252.25, 252.25, 273.59, 273.59)),
  list(c(17.80, 23.81, 23.81, 17.80)), list(c(364.52, 364.52, 367.39, 367.39)),
  list(c(17.40, 24.76, 24.76, 17.40)), list(c(349.69, 349.69, 355.25, 355.25)),
  list(c(77.45, 99.91, 99.91, 77.45)), list(c(508.02, 508.02, 516.11, 516.11)),
  list(c(2.36, 8.36, 8.36, 2.36)), list(c(420.74, 420.74, 442.01, 442.01)),
  list(c(228.59, 232.14, 232.14, 228.59)), list(c(194.65, 194.65, 197.92, 197.92)),
  list(c(239.90, 249.32, 249.32, 239.90)), list(c(204.81, 204.81, 212.25, 212.25)),
  list(c(289.86, 308.59, 308.59, 289.86)), list(c(338.29, 338.29, 344.02, 344.02)),
  list(c(213.34, 219.12, 219.12, 213.34)), list(c(253.07, 253.07, 272.49, 272.49)),
  list(c(3.27, 4.82, 4.82, 3.27)), list(c(81.91, 81.91, 101.60, 101.60)),
  list(c(4.09, 5.82, 5.82, 4.09)), list(c(103.05, 103.05, 127.46, 127.46)),
  list(c(228.53, 232.08, 232.08, 228.53)), list(c(16.00, 16.00, 19.68, 19.68)),
  list(c(291.10, 310.60, 310.60, 291.10)), list(c(162.09, 162.09, 165.60, 165.60)),
  list(c(213.41, 218.17, 218.17, 213.41)), list(c(74.44, 74.44, 95.14, 95.14)),
  list(c(19.34, 20.48, 20.48, 19.34)), list(c(194.71, 194.71, 197.46, 197.46))
)

fig1_mask_dots <- tribble(
  ~x,      ~y,      ~r,    ~fill,
   316.99,  142.20,  0.70, "#FEFFF9",
   334.72,  123.00,  0.33, "#FEFFF9"
)

fig1_arrows <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   297.28,  447.81,  249.11,  439.27, 0.7497,
   324.17,  420.43,  256.13,  394.35, 0.7497,
   337.43,  422.52,  317.31,  375.83, 0.7497,
   359.61,  419.32,  355.21,  376.77, 0.7497,
   341.92,  445.69,  354.97,  449.83, 0.7497,
    86.47,  276.28,   38.30,  267.75, 0.7497,
   127.75,  251.30,  107.62,  205.78, 0.7497,
   129.74,  274.07,  143.90,  277.77, 0.7497,
   124.72,  291.74,  135.84,  294.04, 0.7497,
   140.55,  300.35,  138.70,  311.58, 0.7497,
   113.58,  250.16,   44.69,  225.04, 0.7497,
    85.57,  447.98,   37.40,  439.44, 0.7497,
   113.83,  421.72,   45.79,  395.65, 0.7497,
   139.11,  470.19,  137.26,  481.42, 0.7497,
   119.79,  419.85,  131.17,  423.14, 0.7497,
   295.98,  273.52,  248.56,  266.69, 0.7497,
   299.19,  217.11,  289.18,  204.94, 0.7497,
   249.97,  295.13,  243.49,  317.84, 0.5561,
   335.91,  114.67,  349.62,  117.18, 0.5561,
   325.55,   73.82,  251.58,   44.46, 0.5561,
   298.68,  100.50,  250.00,   86.67, 0.5561,
   338.82,   74.91,  313.72,   25.87, 0.5561,
   362.17,   73.14,  357.35,   28.90, 0.5561,
   352.28,  122.91,  351.36,  135.64, 0.5561,
   250.89,  118.23,  244.24,  141.53, 0.5561,
   300.80,   39.53,  277.35,   29.29, 0.5561,
   342.70,   93.39,  356.73,   99.27, 0.5561
)

fig1_leaders <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,
    63.36,  100.50,   63.39,  103.91,
    53.02,   99.56,   48.20,   95.90,
   142.86,  112.59,  141.11,  118.36,
    96.22,  100.90,   98.14,  104.58,
    87.84,  104.01,   75.82,  118.31,
   105.57,   37.03,  109.20,   38.92,
   114.11,   50.98,  120.06,   55.21,
   128.79,   63.54,  131.50,   67.33,
   139.20,  122.93,  133.93,  127.75
)

fig1_axes <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   372.90,  155.22,  227.94,  155.28, 1.29,
   376.04,   17.41,  376.09,  151.37, 1.19
)

fig1_silhouettes <- tribble(
  ~species,                     ~x0,    ~y0,    ~x1,    ~y1,    ~flip,
  "Anguilla anguilla", 238.53, 390.96, 254.14, 399.10, "none",
  "Atractosteus spatula", 338.56, 366.80, 372.24, 375.73, "h",
  "Wallago attu", 294.41, 364.04, 328.07, 373.34, "none",
  "Oreochromis andersonii", 236.05, 435.89, 247.38, 442.05, "none",
  "Salmo trutta", 356.20, 446.23, 373.01, 452.66, "none",
  "Huso huso", 123.19, 313.62, 156.87, 321.01, "none",
  "Anguilla anguilla", 27.07, 222.31, 42.68, 230.45, "none",
  "Wallago attu", 86.88, 193.89, 120.53, 203.18, "none",
  "Luciobarbus brachycephalus", 141.58, 290.43, 158.43, 296.87, "none",
  "Oreochromis andersonii", 25.10, 264.48, 36.44, 270.64, "none",
  "Salmo trutta", 145.05, 275.34, 161.85, 281.77, "none",
  "Huso huso", 122.71, 483.21, 156.39, 490.60, "none",
  "Anguilla anguilla", 28.28, 391.59, 43.89, 399.73, "none",
  "Oreochromis andersonii", 24.50, 436.39, 35.84, 442.54, "none",
  "Salmo trutta", 132.07, 420.94, 147.35, 426.79, "none",
  "Dermogenys pusilla", 274.04, 200.65, 292.21, 204.06, "none",
  "Oreochromis andersonii", 235.43, 263.53, 246.77, 269.68, "none",
  "Hemiancistrus medians", 238.31, 316.66, 243.24, 328.23, "none",
  "Huso huso", 337.53, 136.76, 372.09, 144.35, "none",
  "Oreochromis andersonii", 235.09, 83.64, 245.56, 89.32, "none",
  "Anguilla anguilla", 234.39, 41.06, 250.41, 49.41, "none",
  "Dermogenys pusilla", 261.88, 24.39, 280.52, 27.89, "none",
  "Atractosteus spatula", 339.69, 19.41, 374.26, 28.58, "h",
  "Wallago attu", 294.77, 14.83, 329.30, 24.37, "none",
  "Luciobarbus brachycephalus", 353.12, 113.06, 370.42, 119.66, "none",
  "Hemiancistrus medians", 238.02, 139.80, 243.08, 151.68, "none",
  "Salmo trutta", 357.67, 97.08, 374.92, 103.68, "none"
)

fig1_texts <- tribble(
  ~label, ~x, ~y, ~size, ~face, ~col, ~rot,
  "f Game fish", 224.15, 359.82, 7.01, "bold", "black", 0,
  "n = 389", 227.01, 368.02, 5.01, "plain", "black", 0,
  "FRic = 11.40%", 227.01, 374.28, 5.01, "plain", "black", 0,
  "SES = -0.10", 227.01, 380.54, 5.01, "plain", "black", 0,
  "P = 0.458", 227.01, 386.80, 5.01, "plain", "black", 0,
  "PC1 (22.7%)", 286.53, 514.02, 5.01, "plain", "black", 0,
  "PC2 (21.1%)", 216.68, 458.84, 5.01, "plain", "black", 90,
  "c Fisheries", 13.99, 189.74, 7.01, "bold", "black", 0,
  "n = 1,688", 16.46, 198.05, 5.01, "plain", "black", 0,
  "FRic = 35.81%", 16.46, 204.31, 5.01, "plain", "black", 0,
  "SES = 0.32", 16.46, 210.57, 5.01, "plain", "black", 0,
  "P = 0.618", 16.46, 216.83, 5.01, "plain", "black", 0,
  "PC1 (22.7%)", 73.96, 343.69, 5.01, "plain", "black", 0,
  "PC2 (21.1%)", 6.44, 272.52, 5.01, "plain", "black", 90,
  "n = 237", 16.90, 368.20, 5.01, "plain", "black", 0,
  "FRic = 8.08%", 16.90, 374.46, 5.01, "plain", "black", 0,
  "SES = 2.83", 16.90, 380.73, 5.01, "plain", "black", 0,
  "P = 0.998", 16.90, 386.99, 5.01, "plain", "black", 0,
  "e Aquaculture", 13.65, 359.97, 7.01, "bold", "black", 0,
  "PC1 (22.7%)", 73.96, 514.02, 5.01, "plain", "black", 0,
  "PC2 (21.1%)", 6.44, 458.84, 5.01, "plain", "black", 90,
  "n = 1,505", 226.75, 198.38, 5.01, "plain", "black", 0,
  "FRic = 29.05%", 226.75, 204.64, 5.01, "plain", "black", 0,
  "SES = -7.17", 226.75, 210.90, 5.01, "plain", "black", 0,
  "P = 0.001", 226.75, 217.16, 5.01, "plain", "black", 0,
  "d Aquarium trade", 224.27, 189.74, 7.01, "bold", "black", 0,
  "PC1 (22.7%)", 282.78, 343.69, 5.01, "plain", "black", 0,
  "PC2 (21.1%)", 216.64, 272.52, 5.01, "plain", "black", 90,
  "PC2 (21.1%)", 4.56, 102.98, 4.20, "plain", "black", 90,
  "Body mass", 129.46, 111.60, 3.86, "plain", "grey20", 0,
  "Standard", 115.76, 126.59, 3.86, "plain", "grey20", 0,
  "body length", 115.76, 131.41, 3.86, "plain", "grey20", 0,
  "Body", 92.42, 107.67, 3.86, "plain", "grey20", 0,
  "elongation", 92.42, 112.50, 3.86, "plain", "grey20", 0,
  "Caudal", 58.43, 117.71, 3.86, "plain", "grey20", 0,
  "peduncle throttling", 58.43, 122.54, 3.86, "plain", "grey20", 0,
  "Pectoral fin size", 48.22, 107.13, 3.86, "plain", "grey20", 0,
  "Vertical eye position", 47.37, 80.44, 3.86, "plain", "grey20", 0,
  "Relative eye size", 47.96, 51.59, 3.86, "plain", "grey20", 0,
  "Oral gape", 90.74, 34.98, 3.86, "plain", "grey20", 0,
  "position", 90.74, 39.81, 3.86, "plain", "grey20", 0,
  "           Pectoral", 110.42, 54.65, 3.86, "plain", "grey20", 0,
  "fin vertical position", 110.42, 59.48, 3.86, "plain", "grey20", 0,
  "            Relative", 119.77, 66.86, 3.86, "plain", "grey20", 0,
  "maxillary length", 119.77, 71.69, 3.86, "plain", "grey20", 0,
  "Body lateral shape", 30.38, 95.29, 3.86, "plain", "grey20", 0,
  "PC1 (22.7%)", 78.39, 175.45, 4.20, "plain", "black", 0,
  "a PCA correlation circle", 19.33, 11.89, 7.01, "bold", "black", 0,
  "Large bod", 358.42, 164.89, 5.01, "plain", "grey20", 0,
  "y mass", 380.98, 164.89, 5.01, "plain", "grey20", 0,
  "Large mouth", 358.42, 171.15, 5.01, "plain", "grey20", 0,
  "Small head", 358.42, 177.41, 5.01, "plain", "grey20", 0,
  "Upward-oriented mouth", 380.32, 14.95, 5.01, "plain", "grey20", 0,
  "Large eyes", 380.32, 21.21, 5.01, "plain", "grey20", 0,
  "Lateral pectoral fins", 380.32, 27.47, 5.01, "plain", "grey20", 0,
  "n = 2,334", 227.50, 19.07, 5.01, "plain", "black", 0,
  "FRic = 45.22%", 227.50, 25.33, 5.01, "plain", "black", 0,
  "SES = 1.10", 227.50, 31.59, 5.01, "plain", "black", 0,
  "P = 0.863", 227.50, 37.85, 5.01, "plain", "black", 0,
  "Small body mass", 214.99, 164.94, 5.01, "plain", "grey20", 0,
  "Small mouth", 214.99, 171.20, 5.01, "plain", "grey20", 0,
  "Large head", 214.99, 177.46, 5.01, "plain", "grey20", 0,
  "Downward-oriented mouth", 380.60, 145.80, 5.01, "plain", "grey20", 0,
  "Small eyes", 380.60, 152.06, 5.01, "plain", "grey20", 0,
  "Ventral pectoral fins", 380.60, 158.32, 5.01, "plain", "grey20", 0,
  "PC1 (22.7%)", 282.73, 168.38, 5.01, "plain", "black", 0,
  "PC2 (21.1%)", 216.64, 100.83, 5.01, "plain", "black", 90,
  "b All uses", 224.35, 11.36, 6.60, "bold", "black", 0
)

figS1_masks <- tribble(
  ~x, ~y,
  list(c(227.93, 233.63, 233.63, 227.93)), list(c(194.55, 194.55, 199.06, 199.06)),
  list(c(227.21, 233.21, 233.21, 227.21)), list(c(349.68, 349.68, 354.64, 354.64)),
  list(c(229.21, 230.67, 230.67, 229.21)), list(c(364.46, 364.46, 367.12, 367.12)),
  list(c(280.71, 320.14, 320.14, 280.71)), list(c(164.09, 164.09, 171.80, 171.80)),
  list(c(211.87, 218.75, 218.75, 211.87)), list(c(75.11, 75.11, 97.03, 97.03)),
  list(c(290.53, 310.90, 310.90, 290.53)), list(c(163.24, 163.24, 170.62, 170.62)),
  list(c(228.67, 232.66, 232.66, 228.67)), list(c(15.58, 15.58, 19.26, 19.26)),
  list(c(16.27, 22.00, 22.00, 16.27)), list(c(194.90, 194.90, 200.21, 200.21)),
  list(c(17.78, 23.77, 23.77, 17.78)), list(c(364.52, 364.52, 367.39, 367.39)),
  list(c(17.37, 24.72, 24.72, 17.37)), list(c(349.69, 349.69, 355.25, 355.25))
)

figS1_mask_dots <- tribble(
  ~x,      ~y,      ~r,    ~fill,
   334.95,  123.79,  0.55, "#FEFFF9",
   335.70,   65.81,  0.43, "#EBFF9F"
)

figS1_arrows <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   312.95,  289.43,  300.21,  322.56, 0.5561,
   301.00,  275.57,  280.47,  301.30, 0.5561,
   296.36,  263.20,  255.97,  264.95, 0.5561,
   251.41,  412.83,  239.63,  400.11, 0.5561,
   312.09,  454.03,  299.36,  487.16, 0.5561,
   302.82,  438.84,  264.84,  447.66, 0.5561,
   297.71,  423.83,  293.98,  373.03, 0.5561,
   318.39,  432.57,  354.88,  441.67, 0.5561,
   342.89,   60.58,  355.80,   56.57, 0.5561,
   256.42,   68.30,  241.12,   47.53, 0.5561,
   311.96,  109.12,  298.81,  143.34, 0.5561,
   303.54,   92.53,  264.32,  101.63, 0.5561,
   299.82,   77.09,  295.71,   32.94, 0.5561,
   316.35,   85.95,  352.06,   95.46, 0.5561,
   321.25,   78.68,  346.91,   33.44, 0.5561,
   299.62,   94.82,  278.42,  121.39, 0.5561,
   294.84,   82.03,  253.12,   83.85, 0.5561,
   107.34,  263.98,  141.92,  273.19, 0.5561,
   138.17,  234.20,  151.55,  223.38, 0.5561,
    41.02,  244.07,   28.28,  228.69, 0.5561,
   100.96,  289.77,   88.23,  322.90, 0.5561,
    93.69,  271.06,   55.71,  279.87, 0.5561,
   113.57,  256.62,  138.41,  212.81, 0.5561,
    40.67,  413.23,   29.32,  398.67, 0.5561,
   100.61,  458.94,   87.87,  492.06, 0.5561,
   113.21,  425.78,  138.06,  381.98, 0.5561,
   108.50,  432.75,  143.07,  441.96, 0.5561
)

figS1_leaders <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,
   119.72,   78.56,  116.50,   78.41,
    55.27,  107.36,   50.05,  105.05,
    89.21,  111.02,   83.86,  124.07,
   115.88,  102.06,  114.75,  103.35,
    92.95,  109.08,   85.15,  148.10,
    96.12,  107.87,  100.66,  110.59,
    49.88,  125.71,   47.09,  125.07,
   105.38,  142.78,  110.05,  140.79,
    55.09,   50.81,   52.43,   56.16,
    72.53,  113.70,   66.40,  127.03,
    77.17,   91.44,   78.97,   93.30
)

figS1_axes <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   373.32,  156.20,  227.68,  156.26, 1.29,
   376.47,   17.52,  376.53,  152.33, 1.19
)

figS1_silhouettes <- tribble(
  ~species,                     ~x0,    ~y0,    ~x1,    ~y1,    ~flip,
  "Dermogenys pusilla", 236.37, 262.91, 254.51, 266.32, "none",
  "Oreochromis andersonii", 293.16, 323.36, 304.48, 329.52, "none",
  "Hemiancistrus medians", 273.91, 299.55, 278.84, 311.13, "none",
  "Wallago attu", 228.75, 443.47, 262.35, 452.77, "none",
  "Anguilla anguilla", 225.55, 390.09, 241.13, 398.23, "none",
  "Oreochromis andersonii", 292.36, 488.14, 303.68, 494.30, "none",
  "Atractosteus spatula", 282.17, 363.18, 315.80, 372.11, "h",
  "Salmo trutta", 355.49, 438.71, 372.27, 445.14, "none",
  "Wallago attu", 228.05, 97.29, 262.75, 106.89, "none",
  "Dermogenys pusilla", 232.73, 81.69, 251.47, 85.21, "none",
  "Anguilla anguilla", 224.84, 41.17, 240.93, 49.58, "none",
  "Oreochromis andersonii", 292.18, 144.67, 303.87, 151.02, "none",
  "Hemiancistrus medians", 272.69, 120.94, 277.78, 132.90, "none",
  "Huso huso", 335.00, 25.71, 369.73, 33.34, "none",
  "Atractosteus spatula", 282.68, 21.71, 317.41, 30.94, "h",
  "Luciobarbus brachycephalus", 356.88, 51.89, 374.26, 58.54, "none",
  "Salmo trutta", 353.29, 92.76, 370.62, 99.40, "none",
  "Wallago attu", 20.11, 275.57, 53.71, 284.86, "none",
  "Anguilla anguilla", 15.60, 219.72, 31.19, 227.86, "none",
  "Oreochromis andersonii", 81.14, 323.65, 92.46, 329.80, "none",
  "Huso huso", 125.78, 205.02, 159.40, 212.41, "none",
  "Luciobarbus brachycephalus", 145.41, 216.40, 162.23, 222.83, "none",
  "Salmo trutta", 144.27, 270.30, 161.05, 276.73, "none",
  "Anguilla anguilla", 16.18, 389.59, 31.76, 397.73, "none",
  "Oreochromis andersonii", 80.79, 492.91, 92.11, 499.07, "none",
  "Huso huso", 125.40, 374.45, 159.02, 381.84, "none",
  "Salmo trutta", 144.17, 438.63, 160.95, 445.06, "none"
)

figS1_texts <- tribble(
  ~label, ~x, ~y, ~size, ~face, ~col, ~rot,
  "n = 1,505", 226.38, 198.38, 5.00, "plain", "black", 0,
  "FRic = 29.05%", 226.38, 204.64, 5.00, "plain", "black", 0,
  "SES = -7.17", 226.38, 210.90, 5.00, "plain", "black", 0,
  "P = 0.001", 226.38, 217.16, 5.00, "plain", "black", 0,
  "d Aquarium trade", 223.91, 189.74, 7.01, "bold", "black", 0,
  "f Game fish", 223.78, 359.82, 7.01, "bold", "black", 0,
  "n = 389", 226.64, 368.02, 5.00, "plain", "black", 0,
  "FRic = 11.40%", 226.64, 374.28, 5.00, "plain", "black", 0,
  "SES = -0.10", 226.64, 380.54, 5.00, "plain", "black", 0,
  "P = 0.458", 226.64, 386.80, 5.00, "plain", "black", 0,
  "Compress", 358.77, 165.93, 5.00, "plain", "grey20", 0,
  "ed body", 381.56, 165.93, 5.00, "plain", "grey20", 0,
  "Downward-oriented eyes", 358.77, 172.19, 5.00, "plain", "grey20", 0,
  "Small head", 358.77, 178.45, 5.00, "plain", "grey20", 0,
  "Small pectoral fins", 380.77, 15.04, 5.00, "plain", "grey20", 0,
  "Elongated body", 380.77, 21.30, 5.00, "plain", "grey20", 0,
  "n = 2,334", 227.24, 19.19, 5.00, "plain", "black", 0,
  "FRic = 45.22%", 227.24, 25.45, 5.00, "plain", "black", 0,
  "SES = 1.10", 227.24, 31.71, 5.00, "plain", "black", 0,
  "P = 0.863", 227.24, 37.97, 5.00, "plain", "black", 0,
  "Elongated body", 214.67, 165.98, 5.00, "plain", "grey20", 0,
  "Upward-oriented eyes", 214.67, 172.24, 5.00, "plain", "grey20", 0,
  "Large head", 214.67, 178.50, 5.00, "plain", "grey20", 0,
  "Large pectoral fins", 381.05, 146.72, 5.00, "plain", "grey20", 0,
  "Compressed body", 381.05, 152.99, 5.00, "plain", "grey20", 0,
  "PC3 (17%)", 292.06, 168.76, 4.09, "plain", "black", 0,
  "PC4 (13.2%)", 215.19, 95.24, 4.09, "plain", "black", 90,
  "b All uses", 224.08, 11.43, 7.01, "bold", "black", 0,
  "c Fisheries", 13.96, 189.74, 7.01, "bold", "black", 0,
  "n = 1,688", 16.43, 198.05, 5.00, "plain", "black", 0,
  "FRic = 35.81%", 16.43, 204.31, 5.00, "plain", "black", 0,
  "SES = 0.32", 16.43, 210.57, 5.00, "plain", "black", 0,
  "P = 0.618", 16.43, 216.83, 5.00, "plain", "black", 0,
  "n = 237", 16.87, 368.20, 5.00, "plain", "black", 0,
  "FRic = 8.08%", 16.87, 374.46, 5.00, "plain", "black", 0,
  "SES = 2.83", 16.87, 380.73, 5.00, "plain", "black", 0,
  "P = 0.998", 16.87, 386.99, 5.00, "plain", "black", 0,
  "e Aquaculture", 13.63, 359.97, 7.01, "bold", "black", 0,
  "b) All usages", 138.89, 17.03, 5.49, "plain", "black", 0,
  "Body", 75.68, 127.71, 4.12, "plain", "grey20", 0,
  "mass", 75.68, 132.86, 4.12, "plain", "grey20", 0,
  "Standard body length", 67.34, 151.56, 4.12, "plain", "grey20", 0,
  "Body elongation", 53.90, 50.09, 4.12, "plain", "grey20", 0,
  "Caudal", 116.90, 98.70, 4.12, "plain", "grey20", 0,
  "peduncle throttling", 116.90, 103.85, 4.12, "plain", "grey20", 0,
  "Pectoral", 110.89, 138.74, 4.12, "plain", "grey20", 0,
  "fin size", 110.89, 143.89, 4.12, "plain", "grey20", 0,
  "Vertical eye", 37.06, 118.78, 4.12, "plain", "grey20", 0,
  "position", 37.06, 123.93, 4.12, "plain", "grey20", 0,
  "Relative eye", 119.86, 77.26, 4.12, "plain", "grey20", 0,
  "size", 119.86, 82.41, 4.12, "plain", "grey20", 0,
  "Oral gape", 101.29, 110.84, 4.12, "plain", "grey20", 0,
  "   position", 101.29, 115.99, 4.12, "plain", "grey20", 0,
  " Pectoral", 53.93, 132.04, 4.12, "plain", "grey20", 0,
  "fin vertical", 53.93, 137.20, 4.12, "plain", "grey20", 0,
  "position", 53.93, 142.35, 4.12, "plain", "grey20", 0,
  "Relative maxillary length", 31.98, 92.51, 4.12, "plain", "grey20", 0,
  "Body lateral shape", 30.21, 103.61, 4.12, "plain", "grey20", 0,
  "a PCA correlation circle", 7.90, 10.78, 7.01, "bold", "black", 0
)

draw_fs_figure <- function(fig, fs_png, circle_png, silhouettes) {
  layout <- function(name) get(paste0(fig, "_", name))
  texts  <- layout("texts")
  under  <- texts$label == "b) All usages"

  function() {
    fig_texts(texts[under, ])

    for (k in which(fs_images$fig == fig)) {
      im <- fs_images[k, ]
      fig_image(fs_png[[fig]][[im$panel]], im$x, im$y, im$w, im$h)
    }
    cp <- circle_place[[fig]]
    fig_image(circle_png[[fig]], cp$box[1], cp$box[2], cp$box[3], cp$box[4], clip = cp$clip)

    masks <- layout("masks")
    for (k in seq_len(nrow(masks))) fig_line(unlist(masks$x[[k]]), unlist(masks$y[[k]]), col = NA, fill = "white")
    dots <- layout("mask_dots")
    for (k in seq_len(nrow(dots))) fig_circle(dots$x[k], dots$y[k], dots$r[k], col = NA, fill = dots$fill[k])

    leaders <- layout("leaders")
    for (k in seq_len(nrow(leaders))) {
      fig_line(c(leaders$x0[k], leaders$x1[k]), c(leaders$y0[k], leaders$y1[k]), col = "#808080", lwd = 0.1251)
    }

    arrows <- layout("arrows")
    for (k in seq_len(nrow(arrows))) {
      a <- arrows[k, ]
      fig_arrow(a$x0, a$y0, a$x1, a$y1, a$lwd)
    }

    axes <- layout("axes")
    for (k in seq_len(nrow(axes))) {
      a <- axes[k, ]
      fig_arrow(a$x0, a$y0, a$x1, a$y1, a$lwd, col = "#BF4808", dashed = FALSE, dot = FALSE, both_ends = TRUE)
    }

    sil <- layout("silhouettes")
    for (k in seq_len(nrow(sil))) {
      s <- sil[k, ]
      fig_silhouette(silhouettes[[s$species]], s$x0, s$y0, s$x1, s$y1, s$flip)
    }

    fig_texts(texts[!under, ])
  }
}

# --- Figures 2 and S3 ---------------------------------------------------------

deficit_images <- tribble(
  ~fig,    ~catg,          ~x,      ~y,      ~w,      ~h,
  "fig2",  "All uses",     194.917,   0.000, 170.315, 170.180,
  "fig2",  "Fisheries",      0.000, 165.098, 170.315, 170.180,
  "fig2",  "Aquarium",     194.917, 165.100, 170.315, 170.180,
  "fig2",  "Game fish",    194.917, 326.820, 170.315, 170.180,
  "fig2",  "Aquaculture",    0.000, 326.820, 170.315, 170.180,
  "figS3", "All uses",     171.160,   0.000, 170.254, 170.207,
  "figS3", "Fisheries",      0.000, 171.513, 170.254, 170.207,
  "figS3", "Aquarium",     171.160, 171.513, 170.254, 170.207,
  "figS3", "Game fish",    171.160, 342.793, 170.254, 170.207,
  "figS3", "Aquaculture",    0.001, 342.791, 170.254, 170.206
)

row_place <- list(
  fig2  = list(box = c(-19.165, -41.805, 419.051, 214.164), clip = c(376.84, 7.82, 387.81, 148.47)),
  figS3 = list(box = c(-46.227, -41.342, 418.360, 213.920), clip = c(349.12, 8.22, 360.08, 148.72))
)

fig2_masks <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,
   278.90,    2.73,  304.61,    9.78,
    84.42,  169.19,  110.14,  176.24,
   278.79,  168.65,  332.05,  175.69,
   269.57,  330.69,  322.82,  337.73,
    69.85,  329.09,  123.11,  336.14,
   376.38,   10.50,  376.84,   10.62,
   387.69,   13.98,  395.86,  144.53
)

fig2_axes <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   352.59,  143.46,  225.91,  143.52, 1.28,
   221.58,   21.91,  221.63,  139.16, 1.19
)

fig2_texts <- tribble(
  ~label, ~x, ~y, ~size, ~face, ~col,
  "a All uses", 221.14, 10.75, 7.01, "bold", "black",
  "\u00b5", 226.29, 18.61, 5.00, "plain", "black",
  "PC2", 229.17, 19.61, 3.25, "plain", "black",
  " = -0.08", 235.50, 18.61, 5.00, "plain", "black",
  "SES = -3.08", 226.29, 24.86, 5.00, "plain", "black",
  "P = 0.002", 226.29, 31.11, 5.00, "plain", "black",
  "\u00b5", 331.57, 125.60, 5.00, "plain", "black",
  "PC1", 334.45, 126.60, 3.25, "plain", "black",
  " = 0.52", 340.78, 125.60, 5.00, "plain", "black",
  "SES = 17.99", 331.57, 131.85, 5.00, "plain", "black",
  "P = 0.999", 331.57, 138.11, 5.00, "plain", "black",
  "b Fisheries", 26.07, 175.85, 7.01, "bold", "black",
  "\u00b5", 29.25, 183.72, 5.00, "plain", "black",
  "PC2", 32.14, 184.72, 3.25, "plain", "black",
  " = -0.02", 38.47, 183.72, 5.00, "plain", "black",
  "SES = -0.54", 29.25, 189.97, 5.00, "plain", "black",
  "P = 0.304", 29.25, 196.22, 5.00, "plain", "black",
  "\u00b5", 135.12, 293.31, 5.00, "plain", "black",
  "PC1", 138.00, 294.31, 3.25, "plain", "black",
  " = 0.54", 144.33, 293.31, 5.00, "plain", "black",
  "SES = 15.45", 135.12, 299.56, 5.00, "plain", "black",
  "P = 0.999", 135.12, 305.82, 5.00, "plain", "black",
  "c Aquarium trade", 221.12, 174.94, 7.01, "bold", "black",
  "\u00b5", 226.29, 183.72, 5.00, "plain", "black",
  "PC2", 229.17, 184.72, 3.25, "plain", "black",
  " = 0.20", 235.50, 183.72, 5.00, "plain", "black",
  "SES = 5.65", 226.29, 189.98, 5.00, "plain", "black",
  "P = 0.999", 226.29, 196.23, 5.00, "plain", "black",
  "\u00b5", 331.57, 293.31, 5.00, "plain", "black",
  "PC1", 334.45, 294.31, 3.25, "plain", "black",
  " = -0.27", 340.78, 293.31, 5.00, "plain", "black",
  "SES = -7.57", 331.57, 299.56, 5.00, "plain", "black",
  "P = 0.001", 331.57, 305.82, 5.00, "plain", "black",
  "e Game fish", 221.12, 337.40, 7.01, "bold", "black",
  "\u00b5", 226.29, 345.41, 5.00, "plain", "black",
  "PC2", 229.17, 346.41, 3.25, "plain", "black",
  " = -0.96", 235.50, 345.41, 5.00, "plain", "black",
  "SES = -12.30", 226.29, 351.66, 5.00, "plain", "black",
  "P = 0.001", 226.29, 357.92, 5.00, "plain", "black",
  "\u00b5", 331.57, 455.05, 5.00, "plain", "black",
  "PC1", 334.45, 456.05, 3.25, "plain", "black",
  " = 2.50", 340.78, 455.05, 5.00, "plain", "black",
  "SES = 31.75", 331.57, 461.31, 5.00, "plain", "black",
  "P = 0.999", 331.57, 467.56, 5.00, "plain", "black",
  "d Aquaculture", 26.23, 336.66, 7.01, "bold", "black",
  "\u00b5", 29.24, 345.40, 5.00, "plain", "black",
  "PC2", 32.12, 346.40, 3.25, "plain", "black",
  " = -0.91", 38.46, 345.40, 5.00, "plain", "black",
  "SES = -9.53", 29.24, 351.65, 5.00, "plain", "black",
  "P = 0.001", 29.24, 357.90, 5.00, "plain", "black",
  "\u00b5", 135.03, 455.07, 5.00, "plain", "black",
  "PC1", 137.92, 456.07, 3.25, "plain", "black",
  " = 2.47", 144.25, 455.07, 5.00, "plain", "black",
  "SES = 24.89", 135.03, 461.32, 5.00, "plain", "black",
  "P = 0.999", 135.03, 467.57, 5.00, "plain", "black",
  "100%", 388.06, 14.21, 3.98, "plain", "black",
  "75%", 388.06, 46.85, 3.98, "plain", "black",
  "50%", 388.06, 78.93, 3.98, "plain", "black",
  "25%", 388.06, 112.14, 3.98, "plain", "black",
  "0%", 388.33, 144.85, 3.98, "plain", "black",
  "Large body mass", 343.27, 152.49, 5.00, "plain", "grey20",
  "Large mouth", 343.27, 158.74, 5.00, "plain", "grey20",
  "Small head", 343.27, 164.99, 5.00, "plain", "grey20",
  "Upward-oriented mouth", 158.33, 16.53, 5.00, "plain", "grey20",
  "Large eyes", 158.33, 22.78, 5.00, "plain", "grey20",
  "Lateral pectoral fins", 158.33, 29.03, 5.00, "plain", "grey20",
  "Small body mass", 206.26, 152.53, 5.00, "plain", "grey20",
  "Small mouth", 206.26, 158.79, 5.00, "plain", "grey20",
  "Large head", 206.26, 165.04, 5.00, "plain", "grey20",
  "Downward-oriented mouth", 158.67, 138.68, 5.00, "plain", "grey20",
  "Small eyes", 158.67, 144.93, 5.00, "plain", "grey20",
  "Ventral pectoral fins", 158.67, 151.19, 5.00, "plain", "grey20"
)

figS3_masks <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,
   200.78,  149.37,  209.40,  153.28,
   314.57,  149.36,  323.19,  153.27,
   255.11,    2.73,  280.82,    9.78,
    84.39,  175.61,  110.10,  182.66,
   255.16,  174.69,  308.39,  181.75,
   245.78,  346.66,  299.02,  353.71,
    69.83,  345.07,  123.07,  352.12,
   348.67,   10.91,  349.12,   11.02,
   359.96,   14.38,  368.12,  144.79
)

figS3_axes <- tribble(
  ~x0,     ~y0,     ~x1,     ~y1,     ~lwd,
   330.24,  143.86,  202.12,  143.92, 1.28,
   197.73,   20.86,  197.79,  139.51, 1.19
)

figS3_texts <- tribble(
  ~label, ~x, ~y, ~size, ~face, ~col,
  "a All uses", 197.38, 10.75, 7.01, "bold", "black",
  "b Fisheries", 26.06, 182.26, 7.01, "bold", "black",
  "c Aquarium trade", 197.35, 181.35, 7.01, "bold", "black",
  "e Game fish", 197.35, 353.38, 7.01, "bold", "black",
  "d Aquaculture", 26.22, 352.63, 7.01, "bold", "black",
  "100%", 360.33, 14.61, 3.98, "plain", "black",
  "75%", 360.33, 47.21, 3.98, "plain", "black",
  "50%", 360.33, 79.26, 3.98, "plain", "black",
  "25%", 360.33, 112.43, 3.98, "plain", "black",
  "0%", 360.59, 145.10, 3.98, "plain", "black",
  "Compressed body", 320.83, 153.13, 5.00, "plain", "grey20",
  "Downward-oriented eyes", 320.83, 159.39, 5.00, "plain", "grey20",
  "Small head", 320.83, 165.64, 5.00, "plain", "grey20",
  "Small pectoral fins", 151.27, 15.42, 5.00, "plain", "grey20",
  "Elongated body", 151.27, 21.67, 5.00, "plain", "grey20",
  "Elongated body", 182.23, 153.07, 5.00, "plain", "grey20",
  "Upward-oriented eyes", 182.23, 159.33, 5.00, "plain", "grey20",
  "Large head", 182.23, 165.58, 5.00, "plain", "grey20",
  "Large pectoral fins", 151.24, 139.01, 5.00, "plain", "grey20",
  "Compressed body", 151.24, 145.27, 5.00, "plain", "grey20",
  "\u00b5", 201.08, 17.53, 5.00, "plain", "black",
  "PC4", 203.97, 18.53, 3.25, "plain", "black",
  " = -0.30", 210.30, 17.53, 5.00, "plain", "black",
  "SES = -13.77", 201.08, 23.78, 5.00, "plain", "black",
  "P", 201.08, 30.05, 5.00, "plain", "black",
  " =", 204.42, 30.05, 5.00, "italic", "black",
  " 0.001", 208.74, 30.05, 5.00, "plain", "black",
  "\u00b5", 29.99, 188.63, 5.00, "plain", "black",
  "PC4", 32.88, 189.63, 3.25, "plain", "black",
  " = -0.31", 39.21, 188.63, 5.00, "plain", "black",
  "SES = -11.84", 29.99, 194.89, 5.00, "plain", "black",
  "P", 29.99, 201.16, 5.00, "plain", "black",
  "=", 34.63, 201.16, 5.00, "italic", "black",
  " 0.001", 37.56, 201.16, 5.00, "plain", "black",
  "\u00b5", 29.23, 359.87, 5.00, "plain", "black",
  "PC4", 32.11, 360.87, 3.25, "plain", "black",
  " = -0.50", 38.44, 359.87, 5.00, "plain", "black",
  "SES = -6.36", 29.23, 366.12, 5.00, "plain", "black",
  "P", 29.23, 372.39, 5.00, "plain", "black",
  " =", 32.57, 372.39, 5.00, "italic", "black",
  " 0.001", 36.88, 372.39, 5.00, "plain", "black",
  "\u00b5", 200.35, 359.88, 5.00, "plain", "black",
  "PC4", 203.24, 360.88, 3.25, "plain", "black",
  " = -0.55", 209.57, 359.88, 5.00, "plain", "black",
  "SES = -9.28", 200.35, 366.14, 5.00, "plain", "black",
  "P", 200.35, 372.41, 5.00, "plain", "black",
  " =", 203.69, 372.41, 5.00, "italic", "black",
  " 0.001", 208.01, 372.41, 5.00, "plain", "black",
  "\u00b5", 200.39, 188.64, 5.00, "plain", "black",
  "PC4", 203.27, 189.64, 3.25, "plain", "black",
  " = -0.36", 209.60, 188.64, 5.00, "plain", "black",
  "SES = -12.45", 200.39, 194.89, 5.00, "plain", "black",
  "P", 200.39, 201.16, 5.00, "plain", "black",
  " =", 203.72, 201.16, 5.00, "italic", "black",
  " 0.001", 208.04, 201.16, 5.00, "plain", "black",
  "\u00b5", 308.30, 126.97, 5.00, "plain", "black",
  "PC3", 311.18, 127.97, 3.25, "plain", "black",
  " = 0.25", 317.51, 126.97, 5.00, "plain", "black",
  "SES = 10.68", 308.30, 133.22, 5.00, "plain", "black",
  "P", 308.30, 139.49, 5.00, "plain", "black",
  " =", 311.64, 139.49, 5.00, "italic", "black",
  " 0.999", 315.95, 139.49, 5.00, "plain", "black",
  "\u00b5", 310.77, 299.74, 5.00, "plain", "black",
  "PC3", 313.65, 300.75, 3.25, "plain", "black",
  " = 0.29", 319.98, 299.74, 5.00, "plain", "black",
  "SES = 8.88", 310.77, 306.00, 5.00, "plain", "black",
  "P", 310.77, 312.27, 5.00, "plain", "black",
  " =", 314.10, 312.27, 5.00, "italic", "black",
  " 0.999", 318.42, 312.27, 5.00, "plain", "black",
  "\u00b5", 310.76, 471.04, 5.00, "plain", "black",
  "PC3", 313.65, 472.05, 3.25, "plain", "black",
  " = 0.14", 319.97, 471.04, 5.00, "plain", "black",
  "SES = 2.06", 310.76, 477.30, 5.00, "plain", "black",
  "P", 310.76, 483.57, 5.00, "plain", "black",
  " =", 314.10, 483.57, 5.00, "italic", "black",
  " 0.979", 318.41, 483.57, 5.00, "plain", "black",
  "\u00b5", 139.49, 471.06, 5.00, "plain", "black",
  "PC3", 142.37, 472.06, 3.25, "plain", "black",
  " = 0.30", 148.70, 471.06, 5.00, "plain", "black",
  "SES = 3.54", 139.49, 477.31, 5.00, "plain", "black",
  "P", 139.49, 483.58, 5.00, "plain", "black",
  " =", 142.83, 483.58, 5.00, "italic", "black",
  " 0.999", 147.14, 483.58, 5.00, "plain", "black",
  "\u00b5", 140.32, 299.74, 5.00, "plain", "black",
  "PC3", 143.21, 300.75, 3.25, "plain", "black",
  " = 0.26", 149.54, 299.74, 5.00, "plain", "black",
  "SES = 8.70", 140.32, 306.00, 5.00, "plain", "black",
  "P", 140.32, 312.27, 5.00, "plain", "black",
  " =", 143.66, 312.27, 5.00, "italic", "black",
  " 0.999", 147.98, 312.27, 5.00, "plain", "black"
)

fig2_bar <- function() {
  grid::grid.rect(
    x = grid::unit(377.242, "native"), y = grid::unit(12.039, "native"),
    width = grid::unit(9.766, "bigpts"), height = grid::unit(131.922, "bigpts"), just = c("left", "top"),
    gp = grid::gpar(col = NA, fill = grid::linearGradient(
      c("#0E1107", "#F1DCB4", "#A44100"), stops = c(0, 0.503676, 1),
      x1 = 0.5, y1 = 0.00069, x2 = 0.5, y2 = 0.99785
    ))
  )
  for (y0 in c(12.039, 44.83, 77.50, 110.35, 143.07)) {
    fig_line(c(387.00, 387.812, 387.812, 387.00), c(y0, y0, y0 + 0.891, y0 + 0.891), col = NA, fill = "#0E1107")
  }
}

draw_deficit_figure <- function(fig, panel_files, row_file) {
  layout <- function(name) get(paste0(fig, "_", name))
  function() {
    for (k in which(deficit_images$fig == fig)) {
      im <- deficit_images[k, ]
      fig_image(panel_files[[fig]][[im$catg]], im$x, im$y, im$w, im$h)
    }
    masks <- layout("masks")
    for (k in seq_len(nrow(masks))) {
      m <- masks[k, ]
      fig_line(c(m$x0, m$x1, m$x1, m$x0), c(m$y0, m$y0, m$y1, m$y1), col = NA, fill = "white")
    }
    rp <- row_place[[fig]]
    fig_image(row_file, rp$box[1], rp$box[2], rp$box[3], rp$box[4], clip = rp$clip, flip_y = TRUE)
    if (fig == "fig2") fig2_bar()

    axes <- layout("axes")
    for (k in seq_len(nrow(axes))) {
      a <- axes[k, ]
      fig_arrow(a$x0, a$y0, a$x1, a$y1, a$lwd, col = "#561531", dashed = FALSE, dot = FALSE, both_ends = TRUE)
    }
    fig_texts(layout("texts"))
  }
}

# --- Figure 3 -----------------------------------------------------------------

fig3_label_pos <- list(
  "All uses"    = c(x = 0.476, y = 0.877),
  "Fisheries"   = c(x = 0.036, y = 0.762),
  "Aquarium"    = c(x = 0.025, y = 0.762),
  "Aquaculture" = c(x = 0.036, y = 0.771),
  "Game fish"   = c(x = 0.037, y = 0.771)
)

# --- Figure 4 -----------------------------------------------------------------

fig4_uses <- c(a = "All uses", b = "Fisheries", c = "Aquarium", d = "Aquaculture", e = "Game fish")

fig4_px <- function(x, y) c(x * 510.86 / 3823, y * 225 / 1682)
fig4_box <- function(x, y, w, h) c(fig4_px(x, y), fig4_px(w, h))

fric_place <- tribble(
  ~panel, ~sx,     ~sy,     ~ox,     ~oy,
  "a",    0.98372, 1.11295,    0.81,   8.88,
  "b",    0.49218, 0.49028, 1745.84,  87.22,
  "c",    0.49186, 0.49028, 2616.41,  87.22,
  "d",    0.49218, 0.49069, 1757.84, 833.66,
  "e",    0.49186, 0.49069, 2616.41, 833.66
)

map_place <- tribble(
  ~panel, ~s,     ~ox,     ~oy,
  "a",    0.8666, -335.8,  362.7,
  "b",    0.4129, 1572.1,  217.9,
  "c",    0.4127, 2443.4,  217.5,
  "d",    0.4128, 1582.2,  963.7,
  "e",    0.4126, 2443.2,  963.5
)

delta_place <- tribble(
  ~panel, ~x,   ~y,
  "a",    1366,  603,
  "b",    2320,  734,
  "c",    3192,  734,
  "d",    2330, 1481,
  "e",    3192, 1479
)

tick_labels <- tibble(
  label = c("Current", "-CR", "-EN", "-VU", "-NT", "-DD"),
  x     = c(353.5, 588.5, 838.5, 1088.5, 1336, 1587.5)
)

nt_circles <- tribble(
  ~x,     ~y,     ~r,   ~lwd,
  1345.7,  975.1, 44.2, 5,
  2419.0,  471.7, 20.0, 3,
  3289.0,  348.1, 19.9, 3,
  2430.5, 1043.0, 19.9, 3,
  3287.7, 1082.3, 19.9, 3
)

fig4_red <- "#B22222"

draw_fig4 <- function(fric_png, map_png, legend_png, delta_nt) {
  for (k in seq_len(nrow(fric_place))) {
    f <- fric_place[k, ]
    b <- fig4_box(f$ox, f$oy, f$sx * 1800, f$sy * 1500)
    fig_image(fric_png[[f$panel]], b[1], b[2], b[3], b[4])
  }
  for (k in seq_len(nrow(map_place))) {
    m <- map_place[k, ]
    b <- fig4_box(m$ox, m$oy, m$s * 2000, m$s * 1600)
    fig_image(map_png[[m$panel]], b[1], b[2], b[3], b[4])
  }
  b <- fig4_box(3245.3, -5.0, 1.1065 * 500, 1.1065 * 1600)
  fig_image(legend_png, b[1], b[2], b[3], b[4])

  p1 <- fig4_px(1345.7, 207); p2 <- fig4_px(1345.7, 931)
  fig_line(c(p1[1], p2[1]), c(p1[2], p2[2]), col = fig4_red, lwd = 3 * 0.1337)
  fig_circle(p1[1], p1[2], 9 * 0.1337, col = NA, fill = fig4_red)
  br <- fig4_px(c(536, 536, 1386, 1386), c(1557, 1593.5, 1593.5, 1560))
  fig_line(br[1:4], br[5:8], col = fig4_red, lwd = 6 * 0.1337)
  for (k in seq_len(nrow(nt_circles))) {
    cc <- fig4_px(nt_circles$x[k], nt_circles$y[k])
    fig_circle(cc[1], cc[2], nt_circles$r[k] * 0.1337, col = fig4_red, lwd = nt_circles$lwd[k] * 0.1337)
  }

  for (k in seq_len(nrow(tick_labels))) {
    p <- fig4_px(tick_labels$x[k], 1571)
    fig_text(tick_labels$label[k], p[1], p[2], size = 6.73, hjust = 0.5)
  }

  for (k in seq_len(nrow(delta_place))) {
    d <- delta_place[k, ]
    p <- fig4_px(d$x, d$y)
    fig_text_runs(
      list(
        list("ΔFRic", 5.05),
        list("NT", 3.08, dy = 0.86),
        list(delta_nt$label[delta_nt$Usage == fig4_uses[[d$panel]]], 5.05)
      ),
      p[1], p[2]
    )
  }

  fig_text("b Fisheries",   246.19,  16.58, size = 7, face = "bold")
  fig_text("a All uses",     26.41,  16.58, size = 7, face = "bold")
  fig_text("c Aquarium",    362.45,  16.72, size = 7, face = "bold")
  fig_text("d Aquaculture", 247.93, 116.64, size = 7, face = "bold")
  fig_text("e Game fish",   362.72, 116.24, size = 7, face = "bold")
  fig_text("Quantile change", 466.00, 16.26, size = 6.08)
}

# --- Figure S5 ----------------------------------------------------------------

figS5_box <- c(x = 1.56, y = 0.58, width = 837.45, height = 586.85)

draw_figS5 <- function(plot, label_r2, label_p) {
  fig_plot(plot, figS5_box[["x"]], figS5_box[["y"]], figS5_box[["width"]], figS5_box[["height"]])
  fig_text(label_r2, 641.78, 46.17, size = 12)
  fig_text("p", 701.48, 46.17, size = 12, face = "italic")
  fig_text(label_p, 711.49, 46.17, size = 12)
}
