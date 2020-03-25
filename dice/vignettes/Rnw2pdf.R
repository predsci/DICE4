rm(list=ls())

require(utils)
require(tools)

Sweave("dice.Rnw")
texi2pdf("dice.tex")
command_to_clean='rm dice.aux dice.bbl dice.blg dice.log dice.out dice.toc'
#system(command_to_clean)
command_to_open = "open dice.pdf"
system(command_to_open)
