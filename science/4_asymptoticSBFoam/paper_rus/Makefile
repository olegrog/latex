FILE=paper
SOURCES=$(FILE).tex pics/ info.txt copyright_transfer.doc

all:		paper.zip paper.pdf

paper.zip: 	$(SOURCES)
	zip -r $@ $(SOURCES)

paper.pdf: 	$(FILE).tex pics/*
	latex $(FILE) && dvipdf $(FILE).dvi

