# Pizza.py script to create a box of kinase particles

c = cdata()
c.seed(327874)
c.box("box",0,0,0,10,10,10)
         
c.part("active-receptor",1000,"box")
c.part("inactive-kin1",1000,"box")
c.part("inactive-kin2",1000,"box")
c.part("inactive-kin3",1000,"box")
c.part("ptase-x",1000,"box")

c.unselect("box")
c.write("data.kinase.box")
