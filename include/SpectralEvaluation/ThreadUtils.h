#include <mutex>
#include <list>

namespace novac
{

// GuardedValue represents a value which can be incremented or decremented in a 
// thread safe manner (using mutex)
struct GuardedValue
{
public:
    GuardedValue()
        : m_value(0)
    {}

    int GetValue() const
    {
        return m_value;
    }

    void Zero()
    {
        std::lock_guard<std::mutex> lock(guard);
        m_value = 0;
    }

    void IncrementValue()
    {
        std::lock_guard<std::mutex> lock(guard);
        ++m_value;
    }

    void DecrementValue()
    {
        std::lock_guard<std::mutex> lock(guard);
        --m_value;
    }

private:
    int m_value = 0;
    std::mutex guard;
};

// GuardedList represents a list with thread safe operations (using mutex)
template<class T>
struct GuardedList
{
public:
    GuardedList()
        : m_items()
    {}

    void Clear()
    {
        std::lock_guard<std::mutex> lock(guard);
        m_items.clear();
    }

    void AddItem(T item)
    {
        std::lock_guard<std::mutex> lock(guard);
        m_items.push_back(item);
    }

    /** Attempts to get the first item in the list.
        @return true if the list contained any items and the first item was retrieved, otherwise false. */
    bool PopFront(T& item)
    {
        std::lock_guard<std::mutex> lock(guard);
        if (m_items.size() == 0)
        {
            return false;
        }
        item = m_items.front();
        m_items.pop_front();
        return true;
    }

    void CopyTo(std::vector<T>& dst)
    {
        std::lock_guard<std::mutex> lock(guard);
        dst.clear();
        dst.reserve(m_items.size());
        for (T item : m_items)
        {
            dst.push_back(item);
        }
    }

    size_t Size()
    {
        std::lock_guard<std::mutex> lock(guard);
        return m_items.size();
    }

private:
    std::list<T> m_items;
    std::mutex guard;
};

}  // namespace novac
