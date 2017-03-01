#ifndef BASE_TYEP_H
#define BASE_TYEP_H

struct Para {
    double alpha;
    double beta;

    double backtrack_link_search_alpha;
    double backtrack_link_search_beta;

    Para()
        : alpha(1),
          beta(1.1),
          backtrack_link_search_alpha(0.2),
          backtrack_link_search_beta(0.5)() {}
};

#endif
