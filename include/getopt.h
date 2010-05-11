#ifdef __cplusplus
extern "C" {
#endif

	extern char *optarg;
int getopt(int argc, char **argv, char *opts);


extern int optind;


extern int opterr;

extern int optopt;

#ifdef __cplusplus
}
#endif
