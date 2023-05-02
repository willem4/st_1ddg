OBJECTS = st_1ddg.o node.o state.o element.o enumtypes.o readtest.o readmesh.o initcon.o testfcn.o output.o eulerforw.o rk3.o cranknicolsonwp.o fluxproc.o flux2proc.o flux.o flux2.o physpara.o sourceproc.o source.o updatenodes.o exactriemann.o	

st_1ddg: $(OBJECTS)
	c++ $(OBJECTS) -lm  -o st_1ddg

clean:
	rm $(OBJECTS)

node.o: node.h
state.o: state.h
element.o: element.h
enumtypes.o: enumtypes.h
readtest.o: readtest.h
readmesh.o: readmesh.h
initcon.o: initcon.h
testfcn.o: testfcn.h
output.o: output.h
eulerforw.o: eulerforw.h
rk3.o: rk3.h
cranknicolsonwp.o: cranknicolsonwp.h
fluxproc.o: fluxproc.h
updatenodes.o: updatenodes.h	
flux2proc.o: flux2proc.h
flux.o: flux.h
flux2.o: flux2.h
physpara.o: physpara.h
sourceproc.o: sourceproc.h
source.o: source.h
exactriemann.o: exactriemann.h
