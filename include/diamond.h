#ifndef DIAMOND_H
#define DIAMOND_H

typedef struct {
  float alpha;
  float beta;
  double l_d;
  double r;
} Geometry;

/* typedef struct { */
/*   double p1x; */
/*   double p2x; */
/*   double cx; */
/*   double cy; */
/* } Arc */

double leading_edge(double x, float alpha);

double trailing_edge(double x, float beta, double l_d);

Arc* calculate_arc(double x, const Geometry* geom);




#endif // DIAMOND_H
