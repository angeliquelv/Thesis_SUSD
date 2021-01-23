# Input: rasterfiles of blocks of the case study area with resInd output, to be merged together
# Output: one rasterbrick with resIndSpatial output of case study area + 
# dataframe for statistical analysis with all output per pixel with NA rows removed and aspect 
# and slope appended to it

rm(list = ls())  # clear workspace

# Install packages and load libraries -------------------------------------

# install packages if required
if(!require(rgdal)){install.packages("rgdal")}
if(!require(raster)){install.packages("raster")}

# load libraries
library(rgdal)
library(raster) # package for raster manipulation

# Set working directory ---------------------------------------------------

workdirnoerror <- "G:/Thesis_data/output_resIndSpatial/output_no_RSS_errors"
workdirtry1 <- "G:/Thesis_data/output_resIndSpatial/output_1st_try_RSS_errors"
workdirtry2 <- "G:/Thesis_data/output_resIndSpatial/output_2nd_try_RSS_errors"
workdirtry3 <- "G:/Thesis_data/output_resIndSpatial/output_3th_try_RSS_errors"
workdirstripes <- "G:/Thesis_data/output_resIndSpatial/output_stripes"

# Check output ------------------------------------------------------------

# Somehow the occurence of RSS errors is not stable 
# Sometimes the output contains many RSS errors or even whole bands with NA
# And when running the exact same again, the RSS errors may disappear
# Some pixels seem to always have RSS errors
# When making the input scene smaller, 
# the frequency of RSS errors and especially NA bands seems to diminish
# The outcome of the nr of bpoints seems to be stable. 
# Occassionally there are pixels with one breakpoint more or less in separate runs

# My approach is to run every scenes. Scenes that contain no erros are completed
# I rerun scenes that contain errors. 
# If there is no overlap in the errors / NA bands between the two scenes 
# I merge two scenes to obtain a complete image
# If there is overlap between NA pixels I rerun the scene one third (final) time
# After this the three scenes are merged to obtain the most complete scene possible
# The most complete out of the three is entered first in the merge function

setwd(workdirnoerror)

no_error_files <- list.files(path = workdirnoerror, pattern = '.grd$') 

# kleurenpalette werkt nog niet helemaal, werkt niet voor scenes zonder 0 (zoals 49)
# Test whether the 
no_error_list <- list()
for(i in 1:length(no_error_files)){ 
  output <- brick(no_error_files[i])
  no_error_list[[i]] <- output
  assign(no_error_files[i],output)
  plot(no_error_list[[i]][[1]],colNA="hotpink",main=no_error_files[i],legend=FALSE,col=terrain.colors(n=(maxValue(no_error_list[[i]][[1]])+1),rev=TRUE))
    legend(x="right",legend=0:maxValue(no_error_list[[i]][[1]]),fill=terrain.colors(n=(maxValue(no_error_list[[i]][[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

setwd(workdirtry1)

try_1_files <- list.files(path = workdirtry1, pattern = '.grd$') 
names_try_1 <-paste0(try_1_files,"_1")

try_1_list <- list()
for(i in 1:length(try_1_files)){ 
  output <- brick(try_1_files[i])
  try_1_list[[i]] <- output
  assign(names_try_1[i],output)
  plot(try_1_list[[i]][[1]],colNA="hotpink",main=names_try_1[i],legend=FALSE,col=terrain.colors(n=(maxValue(try_1_list[[i]][[1]])+1),rev=TRUE))
  legend(x="right",legend=0:maxValue(try_1_list[[i]][[1]]),fill=terrain.colors(n=(maxValue(try_1_list[[i]][[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

setwd(workdirtry2)

try_2_files <- list.files(path = workdirtry2, pattern = '.grd$') 
names_try_2 <-paste0(try_2_files,"_2")

try_2_list <- list()
for(i in 1:length(try_2_files)){ 
  output <- brick(try_2_files[i])
  try_2_list[[i]] <- output
  assign(names_try_2[i],output)
  plot(try_2_list[[i]][[1]],colNA="hotpink",main=names_try_2[i],legend=FALSE,col=terrain.colors(n=(maxValue(try_2_list[[i]][[1]])+1),rev=TRUE))
  legend(x="right",legend=0:maxValue(try_2_list[[i]][[1]]),fill=terrain.colors(n=(maxValue(try_2_list[[i]][[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

setwd(workdirtry3)

try_3_files <- list.files(path = workdirtry3, pattern = '.grd$') 
names_try_3 <-paste0(try_3_files,"_3")

try_3_list <- list()
for(i in 1:length(try_3_files)){ 
  output <- brick(try_3_files[i])
  try_3_list[[i]] <- output
  assign(names_try_3[i],output)
  plot(try_3_list[[i]][[1]],colNA="hotpink",main=names_try_3[i],legend=FALSE,col=terrain.colors(n=(maxValue(try_3_list[[i]][[1]])+1),rev=TRUE))
  legend(x="right",legend=0:maxValue(try_3_list[[i]][[1]]),fill=terrain.colors(n=(maxValue(try_3_list[[i]][[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

setwd(workdirtry3)

try_1_list <- list()
names_try_1 <-paste0(try_3_files,"_1")
try_2_list <- list()
names_try_2 <-paste0(try_3_files,"_2")
try_3_list <- list()
names_try_3 <-paste0(try_3_files,"_3")
merged_names <- paste0(try_3_files) 
for(i in 1:length(try_3_files)){ 
  setwd(workdirtry1)
  output <- brick(try_3_files[i])
  try_1_list[[i]] <- output
  assign(names_try_1[i],output)
  setwd(workdirtry2)
  output <- brick(try_3_files[i])
  try_2_list[[i]] <- output
  assign(names_try_2[i],output)
  setwd(workdirtry3)
  output <- brick(try_3_files[i])
  try_3_list[[i]] <- output
  assign(names_try_3[i],output)
  setwd("G:/Thesis_data/output_resIndSpatial/output_merged")
  # those cells that already have a value in the 1st try (BPNumb/rasterlayer 1 contains data)
  # are masked out of the 2nd try. pixels that had a value in either the 1st or 2nd try are masked out of the 3th try
  # the reason being that, when merging NAs are filled per layer: so per layer it is checked whether the pixel contains data 
  # and data from try 1, 2 or 3 is chosen accordingly.
  # However, this should not happen because the data in the rasterlayers is related to each other per pixel
  # So I want to base the decision, which try is selected only on BPNumb. 
  # If a pixel is NA in try 1, the data from try 2 can be used if thats not NA - that is the data from all rasterlayers 
  # example: otherwise it could happen that DBP is 0 (meaning no drought breakpoint), however the indicators that follow
  # which shoul normally be NA, do have a value because that is taken from another try. This is not possible; data shouldnt be mixed
  try_2_list[[i]]<-mask(try_2_list[[i]],try_1_list[[i]][[1]],inverse=TRUE)
  try_3_list[[i]]<-mask(try_2_list[[i]],try_1_list[[i]][[1]],inverse=TRUE)
  try_3_list[[i]]<-mask(try_2_list[[i]],try_2_list[[i]][[1]],inverse=TRUE)
  merged_output<-merge(try_1_list[[i]],try_2_list[[i]],try_3_list[[i]],filename=merged_names[i],overwrite=TRUE)
  assign(merged_names[i],merged_output)
  ## plot
  plot(merged_output[[1]],colNA="hotpink",main=merged_names[i],legend=FALSE,col=terrain.colors(n=(maxValue(merged_output[[1]])+1),rev=TRUE))
  legend(x="right",legend=0:maxValue(merged_output[[1]]),fill=terrain.colors(n=(maxValue(merged_output[[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

# When simply merging the scenes together, rows of NA appear in between rows of blocks
# After much trial and error I found that when first merging the blocks incrementally in columns
# And then merging the columns together, avoids the problem.

for (i in 1:26){
  print(seq(i,650+i,by=26))
}

# [1]   1  27  53  79 105 131 157 183 209 235 261 287 313 339 365 391 417 443 469 495 521 547 573 599 625 651
# [1]   2  28  54  80 106 132 158 184 210 236 262 288 314 340 366 392 418 444 470 496 522 548 574 600 626 652
# [1]   3  29  55  81 107 133 159 185 211 237 263 289 315 341 367 393 419 445 471 497 523 549 575 601 627 653
# [1]   4  30  56  82 108 134 160 186 212 238 264 290 316 342 368 394 420 446 472 498 524 550 576 602 628 654
# [1]   5  31  57  83 109 135 161 187 213 239 265 291 317 343 369 395 421 447 473 499 525 551 577 603 629 655
# [1]   6  32  58  84 110 136 162 188 214 240 266 292 318 344 370 396 422 448 474 500 526 552 578 604 630 656
# [1]   7  33  59  85 111 137 163 189 215 241 267 293 319 345 371 397 423 449 475 501 527 553 579 605 631 657
# [1]   8  34  60  86 112 138 164 190 216 242 268 294 320 346 372 398 424 450 476 502 528 554 580 606 632 658
# [1]   9  35  61  87 113 139 165 191 217 243 269 295 321 347 373 399 425 451 477 503 529 555 581 607 633 659
# [1]  10  36  62  88 114 140 166 192 218 244 270 296 322 348 374 400 426 452 478 504 530 556 582 608 634 660
# [1]  11  37  63  89 115 141 167 193 219 245 271 297 323 349 375 401 427 453 479 505 531 557 583 609 635 661
# [1]  12  38  64  90 116 142 168 194 220 246 272 298 324 350 376 402 428 454 480 506 532 558 584 610 636 662
# [1]  13  39  65  91 117 143 169 195 221 247 273 299 325 351 377 403 429 455 481 507 533 559 585 611 637 663
# [1]  14  40  66  92 118 144 170 196 222 248 274 300 326 352 378 404 430 456 482 508 534 560 586 612 638 664
# [1]  15  41  67  93 119 145 171 197 223 249 275 301 327 353 379 405 431 457 483 509 535 561 587 613 639 665
# [1]  16  42  68  94 120 146 172 198 224 250 276 302 328 354 380 406 432 458 484 510 536 562 588 614 640 666
# [1]  17  43  69  95 121 147 173 199 225 251 277 303 329 355 381 407 433 459 485 511 537 563 589 615 641 667
# [1]  18  44  70  96 122 148 174 200 226 252 278 304 330 356 382 408 434 460 486 512 538 564 590 616 642 668
# [1]  19  45  71  97 123 149 175 201 227 253 279 305 331 357 383 409 435 461 487 513 539 565 591 617 643 669
# [1]  20  46  72  98 124 150 176 202 228 254 280 306 332 358 384 410 436 462 488 514 540 566 592 618 644 670
# [1]  21  47  73  99 125 151 177 203 229 255 281 307 333 359 385 411 437 463 489 515 541 567 593 619 645 671
# [1]  22  48  74 100 126 152 178 204 230 256 282 308 334 360 386 412 438 464 490 516 542 568 594 620 646 672
# [1]  23  49  75 101 127 153 179 205 231 257 283 309 335 361 387 413 439 465 491 517 543 569 595 621 647 673
# [1]  24  50  76 102 128 154 180 206 232 258 284 310 336 362 388 414 440 466 492 518 544 570 596 622 648 674
# [1]  25  51  77 103 129 155 181 207 233 259 285 311 337 363 389 415 441 467 493 519 545 571 597 623 649 675
# [1]  26  52  78 104 130 156 182 208 234 260 286 312 338 364 390 416 442 468 494 520 546 572 598 624 650 676

# When removing the empty scenes the following blocknrs per column remain, merge these 
# [1]                             183 209 235 261                                      
output_col_1<-merge(NDVI_splitted_183_output.grd,NDVI_splitted_209_output.grd,NDVI_splitted_235_output.grd,
                  NDVI_splitted_261_output.grd,filename="output_tinier_col_1.grd")
# [1]                         158 184 210 236 262 288 
output_col_2<-merge(NDVI_splitted_158_output.grd,NDVI_splitted_184_output.grd,NDVI_splitted_210_output.grd,
                  NDVI_splitted_236_output.grd,NDVI_splitted_262_output.grd,NDVI_splitted_288_output.grd,
                  filename="output_tinier_col_2.grd")
# [1]                     133 159 185 211 237 263 289 315 341    
output_col_3<-merge(NDVI_splitted_133_output.grd,NDVI_splitted_159_output.grd,NDVI_splitted_185_output.grd,
                  NDVI_splitted_211_output.grd,NDVI_splitted_237_output.grd,NDVI_splitted_263_output.grd,
                  NDVI_splitted_289_output.grd,NDVI_splitted_315_output.grd,NDVI_splitted_341_output.grd,
                  filename="output_tinier_col_3.grd")
# [1]                     134 160 186 212 238 264 290 316 342    
output_col_4<-merge(NDVI_splitted_134_output.grd,NDVI_splitted_160_output.grd,NDVI_splitted_186_output.grd,
                  NDVI_splitted_212_output.grd,NDVI_splitted_238_output.grd,NDVI_splitted_264_output.grd,
                  NDVI_splitted_290_output.grd,NDVI_splitted_316_output.grd,NDVI_splitted_342_output.grd,
                  filename="output_tinier_col_4.grd")
# [1]                 109 135 161 187 213 239 265 291 317 343 369   
output_col_5<-merge(NDVI_splitted_109_output.grd,NDVI_splitted_135_output.grd,NDVI_splitted_161_output.grd,
                  NDVI_splitted_187_output.grd,NDVI_splitted_213_output.grd,NDVI_splitted_239_output.grd,
                  NDVI_splitted_265_output.grd,NDVI_splitted_291_output.grd,NDVI_splitted_317_output.grd,
                  NDVI_splitted_343_output.grd,NDVI_splitted_369_output.grd,filename="output_tinier_col_5.grd")
# [1]                 110 136 162 188 214 240 266 292 318 344 370 396 
output_col_6<-merge(NDVI_splitted_110_output.grd,NDVI_splitted_136_output.grd,NDVI_splitted_162_output.grd,
                  NDVI_splitted_188_output.grd,NDVI_splitted_214_output.grd,NDVI_splitted_240_output.grd,
                  NDVI_splitted_266_output.grd,NDVI_splitted_292_output.grd,NDVI_splitted_318_output.grd,
                  NDVI_splitted_344_output.grd,NDVI_splitted_370_output.grd,NDVI_splitted_396_output.grd,
                  filename="output_tinier_col_6.grd")
# [1]              85 111 137 163 189 215 241 267 293 319 345 371 397
output_col_7<-merge(NDVI_splitted_85_output.grd,NDVI_splitted_111_output.grd,NDVI_splitted_137_output.grd,
                  NDVI_splitted_163_output.grd,NDVI_splitted_189_output.grd,NDVI_splitted_215_output.grd,
                  NDVI_splitted_241_output.grd,NDVI_splitted_267_output.grd,NDVI_splitted_293_output.grd,
                  NDVI_splitted_319_output.grd,NDVI_splitted_345_output.grd,NDVI_splitted_371_output.grd,
                  NDVI_splitted_397_output.grd,filename="output_tinier_col_7.grd")
# [1]              86 112 138 164 190 216 242 268 294 320 346 372 398 424   
output_col_8<-merge(NDVI_splitted_86_output.grd,NDVI_splitted_112_output.grd,NDVI_splitted_138_output.grd,
                  NDVI_splitted_164_output.grd,NDVI_splitted_190_output.grd,NDVI_splitted_216_output.grd,
                  NDVI_splitted_242_output.grd,NDVI_splitted_268_output.grd,NDVI_splitted_294_output.grd,
                  NDVI_splitted_320_output.grd,NDVI_splitted_346_output.grd,NDVI_splitted_372_output.grd,
                  NDVI_splitted_398_output.grd,NDVI_splitted_424_output.grd,filename="output_tinier_col_8.grd")
# [1]          61  87 113 139 165 191 217 243 269 295 321 347 373 399 425 451 477 503     
output_col_9<-merge(NDVI_splitted_61_output.grd,NDVI_splitted_87_output.grd,NDVI_splitted_113_output.grd,
                  NDVI_splitted_139_output.grd,NDVI_splitted_165_output.grd,NDVI_splitted_191_output.grd,
                  NDVI_splitted_217_output.grd,NDVI_splitted_243_output.grd,NDVI_splitted_269_output.grd,
                  NDVI_splitted_295_output.grd,NDVI_splitted_321_output.grd,NDVI_splitted_347_output.grd,
                  NDVI_splitted_373_output.grd,NDVI_splitted_399_output.grd,NDVI_splitted_425_output.grd,
                  NDVI_splitted_451_output.grd,NDVI_splitted_477_output.grd,NDVI_splitted_503_output.grd,
                  filename="output_tinier_col_9.grd")
# [1]      36  62  88 114 140 166 192 218 244 270 296 322 348 374 400 426 452 478 504 530
output_col_10<-merge(NDVI_splitted_36_output.grd,NDVI_splitted_62_output.grd,NDVI_splitted_88_output.grd,
                  NDVI_splitted_114_output.grd,NDVI_splitted_140_output.grd,NDVI_splitted_166_output.grd,
                  NDVI_splitted_192_output.grd,NDVI_splitted_218_output.grd,NDVI_splitted_244_output.grd,
                  NDVI_splitted_270_output.grd,NDVI_splitted_296_output.grd,NDVI_splitted_322_output.grd,
                  NDVI_splitted_348_output.grd,NDVI_splitted_374_output.grd,NDVI_splitted_400_output.grd,
                  NDVI_splitted_426_output.grd,NDVI_splitted_452_output.grd,NDVI_splitted_478_output.grd,
                  NDVI_splitted_504_output.grd,NDVI_splitted_530_output.grd,filename="output_tinier_col_10.grd")
# [1]  11  37  63  89 115 141 167 193 219 245 271 297 323 349 375 401 427 453 479 505 531        
output_col_11<-merge(NDVI_splitted_11_output.grd,NDVI_splitted_37_output.grd,NDVI_splitted_63_output.grd,
                   NDVI_splitted_89_output.grd,NDVI_splitted_115_output.grd,NDVI_splitted_141_output.grd,
                   NDVI_splitted_167_output.grd,NDVI_splitted_193_output.grd,NDVI_splitted_219_output.grd,
                   NDVI_splitted_245_output.grd,NDVI_splitted_271_output.grd,NDVI_splitted_297_output.grd,
                   NDVI_splitted_323_output.grd,NDVI_splitted_349_output.grd,NDVI_splitted_375_output.grd,
                   NDVI_splitted_401_output.grd,NDVI_splitted_427_output.grd,NDVI_splitted_453_output.grd,
                   NDVI_splitted_479_output.grd,NDVI_splitted_505_output.grd,NDVI_splitted_531_output.grd,
                   filename="output_tinier_col_11.grd")
# [1]  12  38  64  90 116 142 168 194 220 246 272 298 324 350 376 402 428 454 480 506 532 558  
output_col_12<-merge(NDVI_splitted_12_output.grd,NDVI_splitted_38_output.grd,NDVI_splitted_64_output.grd,
                   NDVI_splitted_90_output.grd,NDVI_splitted_116_output.grd,NDVI_splitted_142_output.grd,
                   NDVI_splitted_168_output.grd,NDVI_splitted_194_output.grd,NDVI_splitted_220_output.grd,
                   NDVI_splitted_246_output.grd,NDVI_splitted_272_output.grd,NDVI_splitted_298_output.grd,
                   NDVI_splitted_324_output.grd,NDVI_splitted_350_output.grd,NDVI_splitted_376_output.grd,
                   NDVI_splitted_402_output.grd,NDVI_splitted_428_output.grd,NDVI_splitted_454_output.grd,
                   NDVI_splitted_480_output.grd,NDVI_splitted_506_output.grd,NDVI_splitted_532_output.grd,
                   NDVI_splitted_558_output.grd,filename="output_tinier_col_12.grd")
# [1]  13  39  65  91 117 143 169 195 221 247 273 299 325 351 377 403 429 455 481 507 533 559   
output_col_13<-merge(NDVI_splitted_13_output.grd,NDVI_splitted_39_output.grd,NDVI_splitted_65_output.grd,
                   NDVI_splitted_91_output.grd,NDVI_splitted_117_output.grd,NDVI_splitted_143_output.grd,
                   NDVI_splitted_169_output.grd,NDVI_splitted_195_output.grd,NDVI_splitted_221_output.grd,
                   NDVI_splitted_247_output.grd,NDVI_splitted_273_output.grd,NDVI_splitted_299_output.grd,
                   NDVI_splitted_325_output.grd,NDVI_splitted_351_output.grd,NDVI_splitted_377_output.grd,
                   NDVI_splitted_403_output.grd,NDVI_splitted_429_output.grd,NDVI_splitted_455_output.grd,
                   NDVI_splitted_481_output.grd,NDVI_splitted_507_output.grd,NDVI_splitted_533_output.grd,
                   NDVI_splitted_559_output.grd,filename="output_tinier_col_13.grd")
# [1]  14  40  66  92 118 144 170 196 222 248 274 300 326 352 378 404 430 456 482 508 534 560     612 638 
output_col_14<-merge(NDVI_splitted_14_output.grd,NDVI_splitted_40_output.grd,NDVI_splitted_66_output.grd,
                   NDVI_splitted_92_output.grd,NDVI_splitted_118_output.grd,NDVI_splitted_144_output.grd,
                   NDVI_splitted_170_output.grd,NDVI_splitted_196_output.grd,NDVI_splitted_222_output.grd,
                   NDVI_splitted_248_output.grd,NDVI_splitted_274_output.grd,NDVI_splitted_300_output.grd,
                   NDVI_splitted_326_output.grd,NDVI_splitted_352_output.grd,NDVI_splitted_378_output.grd,
                   NDVI_splitted_404_output.grd,NDVI_splitted_430_output.grd,NDVI_splitted_456_output.grd,
                   NDVI_splitted_482_output.grd,NDVI_splitted_508_output.grd,NDVI_splitted_534_output.grd,
                   NDVI_splitted_560_output.grd,NDVI_splitted_612_output.grd,NDVI_splitted_638_output.grd,
                   filename="output_tinier_col_14.grd")
# [1]  15  41  67  93 119 145 171 197 223 249 275 301 327 353 379 405 431 457 483 509 535 561 587 613 639 665
output_col_15<-merge(NDVI_splitted_15_output.grd,NDVI_splitted_41_output.grd,NDVI_splitted_67_output.grd,
                   NDVI_splitted_93_output.grd,NDVI_splitted_119_output.grd,NDVI_splitted_145_output.grd,
                   NDVI_splitted_171_output.grd,NDVI_splitted_197_output.grd,NDVI_splitted_223_output.grd,
                   NDVI_splitted_249_output.grd,NDVI_splitted_275_output.grd,NDVI_splitted_301_output.grd,
                   NDVI_splitted_327_output.grd,NDVI_splitted_353_output.grd,NDVI_splitted_379_output.grd,
                   NDVI_splitted_405_output.grd,NDVI_splitted_431_output.grd,NDVI_splitted_457_output.grd,
                   NDVI_splitted_483_output.grd,NDVI_splitted_509_output.grd,NDVI_splitted_535_output.grd,
                   NDVI_splitted_561_output.grd,NDVI_splitted_587_output.grd,NDVI_splitted_613_output.grd,
                   NDVI_splitted_639_output.grd,NDVI_splitted_665_output.grd,filename="output_tinier_col_15.grd")
# [1]  16  42  68  94 120 146 172 198 224 250 276 302 328 354 380 406 432 458 484 510 536 562 588 614 640 666
output_col_16<-merge(NDVI_splitted_16_output.grd,NDVI_splitted_42_output.grd,NDVI_splitted_68_output.grd,
                   NDVI_splitted_94_output.grd,NDVI_splitted_120_output.grd,NDVI_splitted_146_output.grd,
                   NDVI_splitted_172_output.grd,NDVI_splitted_198_output.grd,NDVI_splitted_224_output.grd,
                   NDVI_splitted_250_output.grd,NDVI_splitted_276_output.grd,NDVI_splitted_302_output.grd,
                   NDVI_splitted_328_output.grd,NDVI_splitted_354_output.grd,NDVI_splitted_380_output.grd,
                   NDVI_splitted_406_output.grd,NDVI_splitted_432_output.grd,NDVI_splitted_458_output.grd,
                   NDVI_splitted_484_output.grd,NDVI_splitted_510_output.grd,NDVI_splitted_536_output.grd,
                   NDVI_splitted_562_output.grd,NDVI_splitted_588_output.grd,NDVI_splitted_614_output.grd,
                   NDVI_splitted_640_output.grd,NDVI_splitted_666_output.grd,filename="output_tinier_col_16.grd")
# [1]  17  43  69  95 121 147 173 199 225 251 277 303 329 355 381 407 433 459 485 511 537 563 589 615 641 667
output_col_17<-merge(NDVI_splitted_17_output.grd,NDVI_splitted_43_output.grd,NDVI_splitted_69_output.grd,
                   NDVI_splitted_95_output.grd,NDVI_splitted_121_output.grd,NDVI_splitted_147_output.grd,
                   NDVI_splitted_173_output.grd,NDVI_splitted_199_output.grd,NDVI_splitted_225_output.grd,
                   NDVI_splitted_251_output.grd,NDVI_splitted_277_output.grd,NDVI_splitted_303_output.grd,
                   NDVI_splitted_329_output.grd,NDVI_splitted_355_output.grd,NDVI_splitted_381_output.grd,
                   NDVI_splitted_407_output.grd,NDVI_splitted_433_output.grd,NDVI_splitted_459_output.grd,
                   NDVI_splitted_485_output.grd,NDVI_splitted_511_output.grd,NDVI_splitted_537_output.grd,
                   NDVI_splitted_563_output.grd,NDVI_splitted_589_output.grd,NDVI_splitted_615_output.grd,
                   NDVI_splitted_641_output.grd,NDVI_splitted_667_output.grd,filename="output_tinier_col_17.grd")
# [1]  18  44  70  96 122 148 174 200 226 252 278 304 330 356 382 408 434 460 486 512 538          
output_col_18<-merge(NDVI_splitted_18_output.grd,NDVI_splitted_44_output.grd,NDVI_splitted_70_output.grd,
                   NDVI_splitted_96_output.grd,NDVI_splitted_122_output.grd,NDVI_splitted_148_output.grd,
                   NDVI_splitted_174_output.grd,NDVI_splitted_200_output.grd,NDVI_splitted_226_output.grd,
                   NDVI_splitted_252_output.grd,NDVI_splitted_278_output.grd,NDVI_splitted_304_output.grd,
                   NDVI_splitted_330_output.grd,NDVI_splitted_356_output.grd,NDVI_splitted_382_output.grd,
                   NDVI_splitted_408_output.grd,NDVI_splitted_434_output.grd,NDVI_splitted_460_output.grd,
                   NDVI_splitted_486_output.grd,NDVI_splitted_512_output.grd,NDVI_splitted_538_output.grd,
                   filename="output_tinier_col_18.grd")
# [1]  19  45  71  97 123 149 175 201 227 253 279 305     357 383 409 435 461 487 513  
output_col_19<-merge(NDVI_splitted_19_output.grd,NDVI_splitted_45_output.grd,NDVI_splitted_71_output.grd,
                   NDVI_splitted_97_output.grd,NDVI_splitted_123_output.grd,NDVI_splitted_149_output.grd,
                   NDVI_splitted_175_output.grd,NDVI_splitted_201_output.grd,NDVI_splitted_227_output.grd,
                   NDVI_splitted_253_output.grd,NDVI_splitted_279_output.grd,NDVI_splitted_305_output.grd,
                   NDVI_splitted_357_output.grd,NDVI_splitted_383_output.grd,NDVI_splitted_409_output.grd,
                   NDVI_splitted_435_output.grd,NDVI_splitted_461_output.grd,NDVI_splitted_487_output.grd,
                   NDVI_splitted_513_output.grd,filename="output_tinier_col_19.grd")
# [1]  20  46  72  98 124 150 176 202 228     280             384 410 436   
output_col_20<-merge(NDVI_splitted_20_output.grd,NDVI_splitted_46_output.grd,NDVI_splitted_72_output.grd,
                   NDVI_splitted_98_output.grd,NDVI_splitted_124_output.grd,NDVI_splitted_150_output.grd,
                   NDVI_splitted_176_output.grd,NDVI_splitted_202_output.grd,NDVI_splitted_228_output.grd,
                   NDVI_splitted_280_output.grd,NDVI_splitted_384_output.grd,NDVI_splitted_410_output.grd,
                   NDVI_splitted_436_output.grd,filename="output_tinier_col_20.grd")
# [1]  21  47  73  99 125 151 177 203                             411       
output_col_21<-merge(NDVI_splitted_21_output.grd,NDVI_splitted_47_output.grd,NDVI_splitted_73_output.grd,
                   NDVI_splitted_99_output.grd,NDVI_splitted_125_output.grd,NDVI_splitted_151_output.grd,
                   NDVI_splitted_177_output.grd,NDVI_splitted_203_output.grd,NDVI_splitted_411_output.grd,
                   filename="output_tinier_col_21.grd")
# [1]      48  74 100 126 152 178 204       
output_col_22<-merge(NDVI_splitted_48_output.grd,NDVI_splitted_74_output.grd,NDVI_splitted_100_output.grd,
                   NDVI_splitted_126_output.grd,NDVI_splitted_152_output.grd,NDVI_splitted_178_output.grd,
                   NDVI_splitted_204_output.grd,filename="output_tinier_col_22.grd")
# [1]      49  75 101 127 153 179 205    
output_col_23<-merge(NDVI_splitted_49_output.grd,NDVI_splitted_75_output.grd,NDVI_splitted_101_output.grd,
                   NDVI_splitted_127_output.grd,NDVI_splitted_153_output.grd,NDVI_splitted_179_output.grd,
                   NDVI_splitted_205_output.grd,filename="output_tinier_col_23.grd")
# [1]             102 128 154 180 
output_col_24<-merge(NDVI_splitted_102_output.grd,NDVI_splitted_128_output.grd,NDVI_splitted_154_output.grd,
                   NDVI_splitted_180_output.grd,filename="output_tinier_col_24.grd")
# [1]             103 129 155   
output_col_25<-merge(NDVI_splitted_103_output.grd,NDVI_splitted_129_output.grd,NDVI_splitted_155_output.grd,
                   filename="output_tinier_col_25.grd")
# [1]             104 130                                                                                 
output_col_26<-merge(NDVI_splitted_104_output.grd,NDVI_splitted_130_output.grd,filename="output_tinier_col_26.grd")

output_tinier_merged <- merge(output_col_1,output_col_2,output_col_3,output_col_4,output_col_5,output_col_6,output_col_7,
                            output_col_8,output_col_9,output_col_10,output_col_11,output_col_12,output_col_13,output_col_14,
                            output_col_15,output_col_16,output_col_17,output_col_18,output_col_19,output_col_20,output_col_21,
                            output_col_22,output_col_23,output_col_23,output_col_24,output_col_25,output_col_26,filename="output_tinier_merged.grd")

setwd("G:/Thesis_data/output_resIndSpatial/output_64")
files_64 <- list.files(path = "G:/Thesis_data/output_resIndSpatial/output_64", pattern = '.grd$') 
names_files_64 <-paste0(files_64,"64")

files_64_list <- list()
for(i in 1:length(files_64)){ 
  output <- brick(files_64[i])
  files_64_list[[i]] <- output
  assign(names_files_64[i],output)
  plot(files_64_list[[i]][[1]],colNA="hotpink",main=names_files_64[i],legend=FALSE,col=terrain.colors(n=(maxValue(files_64_list[[i]][[1]])+1),rev=TRUE))
  legend(x="right",legend=0:maxValue(files_64_list[[i]][[1]]),fill=terrain.colors(n=(maxValue(files_64_list[[i]][[1]])+1),rev=TRUE),xpd=NA,inset=-0.1)
}

for (i in 1:8){
  print(seq(i,56+i,by=8))
}

# [1]  1  9 17 25 33 41 49 57
# [1]  2 10 18 26 34 42 50 58
# [1]  3 11 19 27 35 43 51 59
# [1]  4 12 20 28 36 44 52 60
# [1]  5 13 21 29 37 45 53 61
# [1]  6 14 22 30 38 46 54 62
# [1]  7 15 23 31 39 47 55 63
# [1]  8 16 24 32 40 48 56 64

# [1]     9 17 25 33 
output_col_1_64 <- merge(NDVI_splitted_9_output.grd64,NDVI_splitted_17_output.grd64,NDVI_splitted_25_output.grd64,
                       NDVI_splitted_33_output.grd64,filename="output_col_1_64.grd")
# [1]    10 18 26 34
output_col_2_64 <- merge(NDVI_splitted_10_output.grd64,NDVI_splitted_18_output.grd64,NDVI_splitted_26_output.grd64,
                       NDVI_splitted_34_output.grd64,filename="output_col_2_64.grd")
# [1]  3 11 19 27 35 43 51
output_col_3_64 <- merge(NDVI_splitted_3_output.grd64,NDVI_splitted_11_output.grd64,NDVI_splitted_19_output.grd64,
                       NDVI_splitted_27_output.grd64,NDVI_splitted_35_output.grd64,NDVI_splitted_43_output.grd64,
                       NDVI_splitted_51_output.grd64,filename="output_col_3_64.grd")
# [1]  4 12 20 28 36 44 52
output_col_4_64 <- merge(NDVI_splitted_4_output.grd64,NDVI_splitted_12_output.grd64,NDVI_splitted_20_output.grd64,
                       NDVI_splitted_28_output.grd64,NDVI_splitted_36_output.grd64,NDVI_splitted_44_output.grd64,
                       NDVI_splitted_52_output.grd64,filename="output_col_4_64.grd")
# [1]  5 13 21 29 37 45 53 61
output_col_5_64 <- merge(NDVI_splitted_5_output.grd64,NDVI_splitted_13_output.grd64,NDVI_splitted_21_output.grd64,
                       NDVI_splitted_29_output.grd64,NDVI_splitted_37_output.grd64,NDVI_splitted_45_output.grd64,
                       NDVI_splitted_53_output.grd64,NDVI_splitted_61_output.grd64,filename="output_col_5_64.grd")
# [1]  6 14 22 30 38 46 54 62
output_col_6_64 <- merge(NDVI_splitted_6_output.grd64,NDVI_splitted_14_output.grd64,NDVI_splitted_22_output.grd64,
                       NDVI_splitted_30_output.grd64,NDVI_splitted_38_output.grd64,NDVI_splitted_46_output.grd64,
                       NDVI_splitted_54_output.grd64,NDVI_splitted_62_output.grd64,filename="output_col_6_64.grd")
# [1]  7 15 23    39 47 
output_col_7_64 <- merge(NDVI_splitted_7_output.grd64,NDVI_splitted_15_output.grd64,NDVI_splitted_23_output.grd64,
                       NDVI_splitted_39_output.grd64,NDVI_splitted_47_output.grd64,filename="output_col_7_64.grd")
# [1]  8 16 24 
output_col_8_64 <- merge(NDVI_splitted_8_output.grd64,NDVI_splitted_16_output.grd64,NDVI_splitted_24_output.grd64,
                       filename="output_col_8_64.grd")

output_merged_64 <- merge(output_col_1_64,output_col_2_64,output_col_3_64,output_col_4_64,output_col_5_64,
                          output_col_6_64,output_col_7_64,output_col_8_64,filename="output_merged_64.grd")

setwd("G:/Thesis_data/output_resIndSpatial")
output_tinier_merged <- mask(output_tinier_merged,output_merged_64[[1]],inverse=TRUE)
output <- merge(output_merged_64,output_tinier_merged)
names(output) <- names(NDVI_splitted_11_output.grd)
writeRaster(output,"output.grd",overwrite=TRUE)

# Create dataframe --------------------------------------------------------

# Load raster with aspect and slope data
aspect_Ounila <- raster("G:/Thesis_data/Statistics/aspect_Ounila.grd")
slope_Ounila <- raster("G:/Thesis_data/Statistics/slope_Ounila.grd")
names(slope_Ounila) <- "slope"

# Reproject aspect and slope data to match resInd output (NVDI)
aspect_Ounila_reprojected <- projectRaster(from=aspect_Ounila,to=output,method='bilinear')
slope_Ounila_reprojected <- projectRaster(from=slope_Ounila,to=output,method='bilinear')

output <- addLayer(output,aspect_Ounila_reprojected,slope_Ounila_reprojected)

# Create data frame from raster object, with xy coordinates
df_output <- as.data.frame(output, xy=TRUE, na.rm=FALSE) # 1725748 pixels

# Only include rows that don't have only NA values (except 2 columns for x and y)
# First remove the rows that are outside the case study area (outside Ounila shapefile), and only have a value for x and y 
df_output <- df_output[rowSums(is.na(df_output)) != (ncol(df_output))-2,] # 813531 pixels (check this with input!!)
# Remove those pixels that had an RSS failure; everything  NA except for x, y, slope and aspect 
df_output <- df_output[rowSums(is.na(df_output)) != (ncol(df_output))-4,]  # 811277 pixels

##Save dataframe
save(df_output,file="G:/Thesis_data/Statistics/df_output.Rda") 

##Create subsets for DBP=1
df_output_dbp <- subset(df_output,DBP==1) # 630661
##Create subsets with MagObsR<=0 for DBP=1 # zo selecteer je alleen breakpoints met een dip
df_testbrick_MagObsN <- subset(df_testbrick_dbp,DBPMagObsR<=0) # 851

save(df_testbrick_MagObsN,file="G:/Thesis_data/Statistics/df_testbrick_MagObsN.Rda") 

