#include "MyThreadPool.h"

// the constructor just launches some amount of workers
MyThreadPool::MyThreadPool(size_t threads)
	: stop(false)
{
	for (size_t i = 0; i < threads; ++i)
		workers.emplace_back(
			[this]
	{
		for (;;)
		{
			std::function<void()> task; // storage for future execution tasks

			{
				std::unique_lock<std::mutex> lock(this->mtx); // Let's create a lock

				// waiting for the event to occur when the queue is not empty
				this->condition.wait(lock,
									 [this] { return this->stop || !this->tasks.empty(); });

				// checking whether the tasks were disassembled and a stop was set
				if (this->stop && this->tasks.empty())
					return;

				task = std::move(this->tasks.front());
				this->tasks.pop();
			}

			task();
		}
	}
	);
}

// the destructor joins all threads
MyThreadPool::~MyThreadPool()
{
	{
		std::unique_lock<std::mutex> lock(mtx);
		stop = true;
	}
	condition.notify_all();
	for (std::thread &worker : workers)
		worker.join();
}