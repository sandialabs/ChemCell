# Pizza.py script to create kinase particles inside a triangulated sphere
# receptors are put on surface of sphere

c = cdata()
c.sphere("region",0,0,0,6.203505)
c.surf("sphere","region")

c.part2d("active-receptor",1000,"sphere")
c.part("inactive-kin1",1000,"sphere")
c.part("inactive-kin2",1000,"sphere")
c.part("inactive-kin3",1000,"sphere")
c.part("ptase-x",1000,"sphere")

c.unselect("region")
c.write("data.kinase.surf")
