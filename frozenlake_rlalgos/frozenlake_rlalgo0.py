import gym
import numpy as np


def train():
    # epsilon-greedy parameter, in terms of percent among the actions for a given state
    eps = .05
    # discount factor applied on updating the value function for a new policy (anything <= 0.5 ensures convergence)
    dfactor = .6667
    # implement value function as lookup table: (state space, action space)
    valfunc = np.random.random((16, 4))
    oldvf = []
    # residuals for tracking policy convergence (track multiple since true convergence may not occur)
    residuals = [6]
    # episode number dictates how much discount is applied to the value on policy update
    episode = 0
    lastwinepi = 0
    # i learned these from observation (just used for informative purposes in output)
    movemap = {0: "left", 1: "down", 2: "right", 3: "up"}
    
    env = gym.make("FrozenLake-v0")
    env.monitor.start("rlalgo0_experiment0", force=True)
    while np.mean(residuals) > .99:
        observation = env.reset()
        env.render()
        steps = 0
        # track all (state, action) pairs so each value can be updated when the episode is over
        moves = []
        while True:
            # apply epsilon-greedy exploration via perturbation of the value function expected returns
            perturbactvals = [v + v * (2 * eps * np.random.random() - eps) for v in valfunc[observation]]
            action = np.random.choice(np.argwhere(perturbactvals >= np.max(perturbactvals)).flatten())
            print("Desired move: " + movemap[action])
            moves.append((observation, action))
            observation, reward, done, info = env.step(action)
            env.render()
            if done:
                print("Episode finished after {} timesteps".format(steps + 1))
                # go through each (unique) move and update the value function
                oldvf = valfunc.flatten()[:]
                for s, a in set(moves):
                    valfunc[s, a] = dfactor * valfunc[s, a] + reward
                if len(residuals) >= 3:
                    del residuals[0]
                if reward != 0:
                    residuals.append(np.linalg.norm(valfunc.flatten() - oldvf))
                    print("# Episodes since last win: " + str(episode - lastwinepi))
                    print("New residual: " + str(residuals[-1]))
                    lastwinepi = episode
                # note that by exiting here, the done states' values will never be updated, but they don't have to be
                episode += 1
                break
            steps += 1
    
    print("Episodes needed for policy convergence: " + str(episode))
    env.monitor.close()
    return valfunc


def validate100(valfunc):
    score = 0
    lastwinepi = 0
    env = gym.make("FrozenLake-v0")
    env.monitor.start("rlalgo0_experiment0", force=True)
    for episode in range(100):
        observation = env.reset()
        env.render()
        steps = 0
        while True:
            action = np.random.choice(np.argwhere(valfunc[observation] >= np.max(valfunc[observation])).flatten())
            observation, reward, done, info = env.step(action)
            env.render()
            if done:
                print("Episode finished after {} timesteps".format(steps + 1))
                if reward != 0:
                    print("# Episodes since last win: " + str(episode - lastwinepi))
                    lastwinepi = episode
                    score += 1
                break
            steps += 1
    
    print("Percent wins in 100 episodes: "+str(score / 100.))
    env.monitor.close()
    return


if __name__ == "__main__":
    valfunc = train()
    print("Begin validation. Press enter to continue")
    input()
    validate100(valfunc)
    