PROJECT=excluded_points

SOURCES=$(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) 
	latex $(PROJECT)
	latex $(PROJECT)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

