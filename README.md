# physics-equations

a class full of mathematical equations used in physics

# i am not a profesional or expert in physics, but i am studying it, and i am on the beginning of my journey, if you found any mistakes, forgive me for my limited knowledge. thanks if you are contributing to the project, adding more mathematical functions will help me of course.

### how does it work:

there is the class PhysicsEquations, all you have to do is just make a new instance of it :
PhysicsEquations p = new PhysicsEquations();
and then you can do this :
p.userInput("add 1 1")

and it will find the function inside the class and call it for you, and return the value of it, or the results, it supports Vector3, double, arrays of Vector3 or double, and more will be added soon.

### how to add a new function:

and working and adding new mathematical function is soo easy, you just have to added it like this:
public double something(){
return 0;
}
you can add parameters too if you need them, for example:
public double something(double a, double b){
return a \* b - a + b;
}
now you can call this function like this:

p.userInput("something 1 1")
and it will return the results : -1
