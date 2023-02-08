yusuf <- read.table(text="
      group ai n1i ci n2i
     Balcon 14  56 15  58
    Clausen 18  66 19  64
Multicentre 15 100 12  95
     Barber 10  52 12  47
     Norris 21 226 24 228
     Kahler  3  38  6  31
    Ledwich  2  20  3  20
       Snow 19  76 15  67
   Fuccella 15 106  9 114
     Briant  5  62  4  57
       Pitt  0   9  0   8
   Lombardo  8 133 11 127
   Thompson  3  48  3  49
     Hutton  0  16  0  13
     Tonkin  1  42  1  46
      Gupta  0  25  3  25
     Barber 14 221 15 228
      Yusuf  0  11  0  11
    Wilcox1  8 259  7 129
    Wilcox2  6 157  4 158
       CPRG  3 177  2 136", header=TRUE)
usethis::use_data(yusuf, overwrite=TRUE)
