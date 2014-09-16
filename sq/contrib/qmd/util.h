#ifndef UTIL_H
#define UTIL_H
int sym2Z(const char *sym, int type);
void Z2sym(int Z, char *sym);
int findf(FILE *fd, int n, char *str, ...);
void Z2kg(struct molecule_t *ptr);
#endif
