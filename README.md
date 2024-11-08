[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/CbeSPAUj)


<center><h1>PROGRAMMING ASSIGNMENT 8</h1>

**Due: Wednesday, November 16, 2024, 2:30 PM**</center>

---
## **Instructions**
1. Click on `Accept Assignment` at the link posted on Moodle. This is equivalent to forking the repository.</font>
2. Navigate to the repository you just forked. It should be named `pa8-<your-github-username>`
3. Click on the green `Code` button. Click on SSH tab and copy the SSH link (it is important to use SSH and not HTTPS).
4. Open the terminal on your computer and navigate to the folder where you want to clone the repository.
5. Run `git clone <copied-ssh-link>`. Replace `<copied-ssh-link>` with the link you copied in step 3.
6. In Visual Studio Code, click on File -> Open and open the folder of the repository. 

## Evaluation

Your submission will be evaluated based on the correctness **_and completeness_** of your implementation and test cases. If your implementation is correct but your test cases are incomplete, you will lose points. If your implementation is incorrect, you will lose even more points.

## **Code Review**
This assignment will only be graded after you schedule and complete an in-person code review 
with the instructor within 5 working days after the due date, if not earlier. \
Use [this calendly link](https://calendly.com/ssultan-dpq/) provided on the course website to schedule a code review.

---

<br/>
<h1><center>Shortest Path and Graph Implementation</h1></center>

<br/>

In this assignment, you are to implement the `Graph` interface in the `Graph.java` file provided in `app/src/main/java/pa8` folder. The interface defines the methods that you need to implement. 

   Make sure to read the java doc comments in the `Graph` interface to understand the requirements of each method.

In total, you are to implement four classes:

1. `GraphMatrix` implements the `Graph` interface 

    * Internally uses a 2D array to represent the adjacency matrix of an _UNWEIGHTED_ graph. 

2. `GraphList` implements the `Graph` interface

    * Internally uses `Map<Integer, ArrayList<Integer>>` to represent the adjacency list of an _UNWEIGHTED_ graph.

    * Use the Java built-in `Map` interface to represent the adjacency list. The key of the map is the vertex number, and the value is an `ArrayList` of adjacent vertices.

    * You can import the `Map` interface using `import java.util.Map;` in your code and the `ArrayList` class using `import java.util.ArrayList;`.

    * You can create a new `Map` object using `Map<Integer, ArrayList<Integer>> map = new HashMap<>();`.

3. `GraphMatrixWeighted` extends `GraphMatrix`

    * Internally uses a 2D array to represent the adjacency matrix of a _WEIGHTED_ graph. 

4. `GraphListWeighted` extends `GraphList` 

    * Internally, uses `Map<Integer, ArrayList<Pair<Integer, Integer>>>` to represent the adjacency list of a _WEIGHTED_ graph. 

    * The `Pair` class is a simple class that stores two values of type `Integer`. You can import the `Pair` class using `import javafx.util.Pair;` in your code.

    * The first value in the pair is the vertex number, and the second value is the weight of the edge between the vertex and the adjacent vertex.

    * You can create a new `Pair` object using `Pair<String, Integer> pair = new Pair<>("exampleKey", 10);`.

    * You can access the first and second values of the pair using `pair.getKey()` and `pair.getValue()`, respectively.

<br/>

Feel free to define any additional <u>private</u> helper methods or classes as needed. 

You are to write JUnit test cases in corresponding test classes for each of the four classes above. Make sure to test all the methods keeping in mind any <u>edge cases</u> that may arise.
