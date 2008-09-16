# Pizza.py script to create a box of A,B particles

c = cdata()
c.seed(573674)
c.box("box",0,0,0,21.54435,21.54435,21.54435)
c.part("A",3000,"box")
c.part("B",1000,"box")
c.write("data.abc.box","A","B")
