#include <memory>
#include <vector>

#define CONNECT(sender, signal, receiver, method) ((sender)->signal.Bind(receiver, method));

template<typename ... Args>
class SlotBase {
public:
    virtual void Exec(Args ...) = 0;
    
};

template<typename T, typename ... Args>
class Slot : public SlotBase<Args ...> {
public:
    Slot(T *obj, void (T::*func)(Args ...)) {
        m_slotTarget = obj;
        m_func = func;
    }
    void Exec(Args ... para) {
        (m_slotTarget->*m_func)(para ...);
    }

private:
    T* m_slotTarget = nullptr;
    void (T::*m_func)(Args ...);
};

template<typename ... Args>
class Signal {
public:
    template<typename T>
    void Bind(T* pObj, void (T::*func)(Args ...)) {
        //auto pp = std::unique_ptr<Slot<T, Args ...>>(new Slot<T, Args ...>(pObj, func));
        //m_slots.push_back(pp);
        m_slots.push_back(new Slot<T, Args ...>(pObj, func));
    }
    void operator()(Args ... para) {
        for (const auto &x : m_slots) {
            x->Exec(para...);
        }
    }

private:
    //std::vector<std::unique_ptr<SlotBase<Args ...>>> m_slots;
    std::vector<SlotBase<Args ...>*> m_slots;
};

