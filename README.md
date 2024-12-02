# Performance of Causal Forests in Identifying Heterogeneous Treatment Effects
---

This is the final project of the Microeconometrics course (Master's) in the summer semester 2021. 

---

The Jupyter Notebook of my project can be accessed here: 

<a href="https://nbviewer.jupyter.org/github/OpenSourceEconomics/ose-data-science-course-projeect-tihaup/blob/master/Forced_Attendance_Notebook.ipynb" 
   target="_parent">
   <img align="center"
  src="https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.png"
      width="109" height="20">
</a>

This project explores the use of causal forests to identify heterogeneous treatment effects. Causal forests, a method inspired by matching estimation, construct counterfactuals by matching observations that fall into the same leaf. The results show that causal forests excel when there is correlation among regressors or when higher-degree polynomials and interactions are added to the data. They also handle selection bias well, making them valuable for observational studies where random treatment assignment is not guaranteed.

Performance tends to decrease as the number of regressors increases, but causal forests maintain robust performance even with large sample sizes. As an application, the Lalonde (1986) dataset was analyzed to uncover treatment effects across different subgroups. It was found that individuals around 30 years old with a higher level of education benefit the most from the NSW job market program, and those with lower income before treatment experience higher treatment effects.

Overall, causal forests prove to be a powerful tool for evaluating the effectiveness of programs tailored to specific subgroups, offering more precise insights than simply estimating the Average Treatment Effect for the entire sample.
