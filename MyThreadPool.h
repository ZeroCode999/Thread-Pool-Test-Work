#ifndef MYTHREADPOOL_H
#define MYTHREADPOOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <atomic>

class MyThreadPool
{
public:
	MyThreadPool(size_t);
	template<class F, class... Args> 	// The signature of this function is std :: async
	auto add_task(F&& f, Args&&... args)
		->std::future<typename std::result_of<F(Args...)>::type>;
	~MyThreadPool();

private:
	// need to keep track of threads so we can join them
	std::vector< std::thread > workers;
	// the task queue
	std::queue< std::function<void()> > tasks;

	// synchronization
	std::mutex mtx;
	std::condition_variable condition;
	
	// atomic is needed so that the value is not cached
	std::atomic<bool> stop;
};

// add new work item to the pool
template<class F, class... Args>
auto MyThreadPool::add_task(F&& f, Args&&... args)
-> std::future<typename std::result_of<F(Args...)>::type>
{
	using return_type = typename std::result_of<F(Args...)>::type;

	/*
	* packaged_task stores tasks with deferred filling, exception handling
	* make_shared is used, for simplicity, because it is copied, this is a mandatory requirement of std :: function for tasks placed in it
	*/
	auto task = std::make_shared< std::packaged_task<return_type()> >(
		std::bind(std::forward<F>(f), std::forward<Args>(args)...)
		);

	std::future<return_type> res = task->get_future();
	{
		std::unique_lock<std::mutex> lock(mtx);

		// don't allow enqueued after stopping the pool
		if (stop)
			throw std::runtime_error("enqueue on stopped ThreadPool");

		/*tasks.emplace([task]() { (*task)(); });*/
		tasks.push([task]() { (*task)(); });
	}
	
	// new task notification
	condition.notify_one();
	return res;
}

#endif

