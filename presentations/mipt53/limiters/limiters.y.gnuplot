set terminal table; set output "limiters.y.table"; set format "%.5f"
set samples 25; plot [x=0:4] (2*y)/(1+y*y)
