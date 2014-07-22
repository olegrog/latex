set terminal table; set output "limiters.Leer.table"; set format "%.5f"
set samples 25; plot [x=0:4] (2*x)/(1+x)
