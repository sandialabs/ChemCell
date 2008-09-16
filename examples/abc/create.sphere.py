# Pizza.py script to create a box of A,B particles

c = cdata()
c.seed(573674)
c.sphere("sphere",0,0,0,13.36505)
c.part("A",3000,"sphere")
c.part("B",1000,"sphere")
c.write("data.abc.sphere","sphere","A","B")
