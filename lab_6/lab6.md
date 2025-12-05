# Lab 6 Functions
Brad Hunter (PID: A69038089)

Example input data

``` r
student1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
student3 <- c(90, NA, NA, NA, NA, NA, NA, NA)
```

Let’s create a function to grade homework

> Q1 Write a function grade() to determine an overall grade from a
> vector of student homework assignment scores dropping the lowest
> single score. If a student misses a homework (i.e. has an NA value)
> this can be used as a score to be potentially dropped. Your final
> function should be adquately explained with code comments and be able
> to work on an example class gradebook such as this one in CSV format:
> “https://tinyurl.com/gradeinput”

``` r
#' Calculates average scores for a vector of homework scores.
#' Drops the lowest single score, NAs are counted as 0.
#'
#' @param d Numeric vector of homework scores
#'
#' @returns Average score
#' @export
#'
#' @examples
#'  student <- c(100, NA, 90, 90, 90, 90, 97, 80)
#'  grade(student)
#'  
grade <- function(d) {
 d[is.na(d)] <- 0
 mean(d[-which.min(d)]) 
}
```

Let’s test our function

``` r
grade(student1)
```

    [1] 100

``` r
grade(student2)
```

    [1] 91

``` r
grade(student3)
```

    [1] 12.85714

Let’s input the test grades of the example grade book and calculate
their grades.

``` r
url <- "https://tinyurl.com/gradeinput"
book <- read.csv(url, row.names=1)
graded <- apply(book,1,grade)
graded
```

     student-1  student-2  student-3  student-4  student-5  student-6  student-7 
         91.75      82.50      84.25      84.25      88.25      89.00      94.00 
     student-8  student-9 student-10 student-11 student-12 student-13 student-14 
         93.75      87.75      79.00      86.00      91.75      92.25      87.75 
    student-15 student-16 student-17 student-18 student-19 student-20 
         78.75      89.50      88.00      94.50      82.75      82.75 

> Q2. Using your grade() function and the supplied gradebook, Who is the
> top scoring student overall in the gradebook?

``` r
names(graded)[which.max(graded)]
```

    [1] "student-18"

The student with the top score is student-18

> Q3. From your analysis of the gradebook, which homework was toughest
> on students (i.e. obtained the lowest scores overall?

``` r
assignments <- apply(book,2,median, na.rm=T)
names(assignments)[which.min(assignments)]
```

    [1] "hw2"

``` r
boxplot(book)
```

![](lab6_files/figure-commonmark/unnamed-chunk-6-1.png)

hw2 was the toughest assignment for the students.

> Q4. Optional Extension: From your analysis of the gradebook, which
> homework was most predictive of overall score (i.e. highest
> correlation with average grade score)?

``` r
book[is.na(book)] <- 0
cor.grade <- apply(book,2,cor, x=graded)
cor.grade
```

          hw1       hw2       hw3       hw4       hw5 
    0.4250204 0.1767780 0.3042561 0.3810884 0.6325982 

The homework that was most predictive of overall score was hw5.
