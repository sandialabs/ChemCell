# Pizza.py script to create kinase particles inside a sphere

c = cdata()
c.sphere("sphere",0,0,0,6.203505)

c.part("active-receptor",1000,"sphere")
c.part("inactive-kin1",1000,"sphere")
c.part("inactive-kin2",1000,"sphere")
c.part("inactive-kin3",1000,"sphere")
c.part("ptase-x",1000,"sphere")

c.unselect("sphere")
c.write("data.kinase.sphere")
