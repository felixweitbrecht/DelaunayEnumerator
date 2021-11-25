JFLAGS = -g
JC = javac

.SUFFIXES: .java .class

.java.class:
	$(JC) $(JFLAGS) $*.java

CLASSES = \
		  delaunayKD/misc/DoSomething.java

default: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) delaunayKD/AllSimplicesFinder.class
	$(RM) delaunayKD/*/*.class
	
run:
	java delaunayKD.misc.DoSomething
