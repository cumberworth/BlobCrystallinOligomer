// movetype.h

#ifndef MOVETYPE_H
#define MOVETYPE_H

namespace movetype {

    class MCMovetype {
        // MC movetype interface
        public:
            virtual bool move() = 0;

    };

    class VMMCMovetype {
        public:
            bool move();
    };
}

#endif // MOVETYPE_H
