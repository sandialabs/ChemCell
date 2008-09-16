# Pizza.py script to create a box of Lotka particles

c = cdata()
c.seed(598945)
c.box("box",0,0,0,200,200,20)
c.part("y1",1000,"box")
c.part("y2",1000,"box")
c.write("data.lotka.box","y1","y2")
